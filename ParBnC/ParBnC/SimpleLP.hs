
import Numeric.IEEE ( IEEE(infinity) )
import Data.List
import Numeric.LinearAlgebra as LA
import Control.Parallel(par, pseq)
import Control.Parallel.Strategies

data MatConstraints = MatVec [[Double]] [Double] deriving Show

data VariableType = INTEGER | CONTINUOUS deriving (Show, Eq)

data ProblemType = LP | MIP deriving Show

data ObjectiveType = Maximization | Minimization deriving (Show, Eq)

epsilonTol :: R
epsilonTol = 1e-6

isInt :: (RealFrac a) => a -> Bool
isInt x = x == fromInteger (round x)

toTableau :: LA.Vector R -> LA.Matrix R -> LA.Vector R -> LA.Matrix R
toTableau costC matA constB = tab where 
    xb = LA.fromColumns $ LA.toColumns matA ++ [constB]
    z = mappend costC $ vector [0]
    tab = LA.fromRows $ LA.toRows xb ++ [z]

costCheck :: ObjectiveType -> (R -> Bool)
costCheck Maximization = (> 0)
costCheck Minimization = (< 0)

boundCheck :: ObjectiveType -> (R -> Bool)
boundCheck Maximization = (< 0)
boundCheck Minimization = (> 0)

isImprovable :: ObjectiveType -> LA.Matrix R -> Bool
isImprovable obj tab = any (costCheck obj) $ LA.toList cost where 
    cost = subVector 0 (cols tab - 1) $ tab ! (rows tab - 1)

isImprovableDual :: ObjectiveType -> LA.Matrix R -> Bool
isImprovableDual obj tab = any (boundCheck obj) $ LA.toList bounds where 
    lastCol = last $ LA.toColumns tab
    bounds = subVector 0 (rows tab - 1) lastCol

getPivotPosition :: ObjectiveType -> LA.Matrix R -> (Int, Int)
getPivotPosition obj tab = (row, column) where 
    z = tab ! (rows tab - 1)
    cost = subVector 0 (cols tab - 1) z
    column = LA.maxIndex cost
    getElem rowEq
        | elem <= 0 = infinity::R
        | otherwise = (rowEq ! (LA.size rowEq - 1)) / elem
        where 
            elem = rowEq ! column
    restrictions = map getElem $ init (LA.toRows tab)

    Just row =  elemIndex (minimum restrictions) restrictions
    
getPivotPositionDual :: ObjectiveType -> LA.Matrix R -> (Int, Int)
getPivotPositionDual obj tab = (row, column) where 
    lastCol = last $ LA.toColumns tab
    bounds = subVector 0 (rows tab - 1) lastCol
    row = head $ LA.find (boundCheck obj) bounds
    getElem rowEq
        | elem >= 0 = infinity::R
        | otherwise = elem / (rowEq ! (LA.size rowEq - 1))
        where 
            elem = rowEq ! row
    restrictions = map getElem $ init (LA.toColumns tab)
    Just column =  elemIndex (minimum restrictions) restrictions

pivotStep :: Matrix R -> (Int, Int) -> Matrix R
pivotStep tab (row, column) = newTableau where
    pivotVal = (tab ! row) ! column
    newPivotRow = cmap (/ pivotVal) $ tab ! row
    updateRow rowIdx rowEq 
        | rowIdx == row = newPivotRow
        | otherwise = rowEq - newPivotRow * LA.scalar (rowEq ! column)
    newTableau = LA.fromRows $ zipWith updateRow [0..] $ LA.toRows tab

isBasic :: Vector R -> Bool
isBasic colVec = (LA.sumElements colVec == 1) && (length zeroVec == colVecLen) where
    zeroVec = filter (== 0) $ LA.toList colVec
    colVecLen = LA.size colVec - 1

getSolution :: Matrix R -> Vector R
getSolution tab = solution where 
    colSize = LA.cols tab
    columns = LA.toColumns $ LA.takeColumns (colSize - 1) tab
    lastCol = LA.flatten $ LA.dropColumns (colSize - 1) tab
    findSol colVec
        | isBasic colVec = sol 
        | otherwise = 0::R where
            oneIndex = head $ LA.find (==1) colVec
            sol = lastCol ! oneIndex
    solution = LA.fromList $ map findSol columns

updateTab :: ObjectiveType -> LA.Matrix R -> LA.Matrix R
updateTab obj tab
    | not $ isImprovable obj tab = tab
    | otherwise = updateTab obj newTab where 
        pivotPos = getPivotPosition obj tab
        newTab = pivotStep tab pivotPos

updateTabDual :: ObjectiveType -> LA.Matrix R -> LA.Matrix R
updateTabDual obj tab 
    | not $ isImprovableDual obj tab = tab
    | otherwise = updateTabDual obj newTab where 
        pivotPos = getPivotPositionDual obj tab
        newTab = pivotStep tab pivotPos

isImprovableMixed :: LA.Matrix R -> Bool
isImprovableMixed tab = any (> 0) $ LA.toList costs where
    lastRow = last $ LA.toRows tab
    costs = subVector 0 (LA.cols tab - 1) lastRow

getPhaseOneTab :: LA.Matrix R -> (LA.Matrix R, LA.Vector R, Bool)
getPhaseOneTab tab = (newTab, oldCost, needed) where 
    (rowSize, colSize) = LA.size tab
    (constRows, [oldCost]) = splitAt (rowSize - 1) $ LA.toRows tab
    updateRow rowVec
        | boundVal < 0 = -rowVec
        | otherwise = rowVec where 
            boundVal = last $ LA.toList rowVec
    zeroRow = LA.konst (0::R) colSize
    mixedRows = filter (\rowVec -> (rowVec ! (LA.size rowVec - 1)) < 0) constRows
    needed = not $ null mixedRows
    sumRow = sum $ map (\x -> -x) $ mixedRows ++ [zeroRow]
    newConstRows = map updateRow constRows
    newTab = LA.fromRows $ newConstRows ++ [sumRow]

getPivotPositionMixed :: LA.Matrix R -> (Int, Int)
getPivotPositionMixed tab = (row, column) where 
    z = tab ! (rows tab - 1)
    colSize = cols tab - 1
    cost = subVector 0 colSize z
    maxCost = LA.maxElement cost
    column = head $ LA.find (== maxCost) cost
    getElem rowEq
        | elem == 0 = infinity::R
        | val < 0 = infinity::R
        | otherwise = val
        where 
            elem = rowEq ! column
            val = (rowEq ! (LA.size rowEq - 1)) / (rowEq ! column)
    restrictions = map getElem $ init (LA.toRows tab)
    minRatio = minimum restrictions
    row =  last [idx | (idx, ratio) <- zip [0..colSize] restrictions, ratio == minRatio]

updateMixedTab :: LA.Matrix R -> LA.Matrix R
updateMixedTab tab
    | not $ isImprovableMixed tab = tab
    | otherwise = updateMixedTab newTab where 
        pivotPos = getPivotPositionMixed tab
        newTab = pivotStep tab pivotPos
    
simplexWithTab :: ObjectiveType -> LA.Matrix R -> (R, LA.Vector R, LA.Matrix R)
simplexWithTab obj tab = (optVal, solution, lastTab) where 
    lastTab = updateTab obj tab
    solution = getSolution lastTab
    (rowSize, colSize) = LA.size lastTab
    lastVal = lastTab ! (rowSize - 1) ! (colSize - 1)
    optVal = if obj == Maximization then lastVal else (-lastVal)

simplexWithMixedTab :: ObjectiveType -> LA.Matrix R -> (R, LA.Vector R, LA.Matrix R)
simplexWithMixedTab obj tab
    | infeasible = (-infinity, getSolution interTab, interTab)
    | otherwise  = (lastVal, solution, lastTab) where
        (phaseOneTab, oldCost, needed) = getPhaseOneTab tab
        interTab
            | needed = updateMixedTab phaseOneTab
            | otherwise = phaseOneTab
        (rowSize, colSize) = LA.size interTab
        lastPhaseOneVal = interTab ! (rowSize - 1) ! (colSize - 1)
        infeasible = not $ isClose lastPhaseOneVal (0::R)

        (constRows, _) = splitAt (rowSize - 1) $ LA.toRows interTab

        phaseTwoTab = LA.fromRows $ constRows ++ [oldCost]
        lastTab = updateTab obj phaseTwoTab
        solution = getSolution lastTab
        lastVal = lastTab ! (rowSize - 1) ! (colSize - 1)
        optVal = lastVal

b :: Vector R
b = LA.fromList [5,28::R]
matA :: Matrix R
matA = LA.fromLists [[1,1,1,0],[4,7,0,1::R]]
c :: Vector R
c = fromList [5,6,0,0::R]

oldTab :: Matrix R
oldTab = toTableau c matA b
tab :: Matrix R
z :: Vector R
val :: R
(val, z, tab) = simplexWithTab Maximization $ toTableau c matA b

newTab :: Matrix R
newTab = addGomoryCut tab $ getGomoryCut $ tab ! 0
getGomoryCut :: Vector R -> Vector R
getGomoryCut rowVec = gomoryCons where
    [varPart, constPart] = takesV [LA.size rowVec - 1, 1] rowVec 
    posDec num = -posFrac where
        (intPart, fracPart) = properFraction num
        posFrac = if fracPart >= 0.0 then fracPart else 1.0 + fracPart
    gomoryCons = vjoin [cmap posDec varPart, vector [1::R], cmap posDec constPart]

addSlackColumn :: Matrix R -> Matrix R
addSlackColumn tab = newTab where
    columns = LA.toColumns tab
    (rowSize, colSize) = LA.size tab
    (varColumns, lastColumn) = splitAt (colSize - 1) columns
    zeroVec = LA.konst 0 rowSize :: Vector R
    extendedTab = varColumns ++ [zeroVec] ++ lastColumn
    newTab = LA.fromColumns extendedTab

addNewRow :: Matrix R -> Vector R -> Matrix R
addNewRow tab newRow = newTab where
    rows = LA.toRows tab
    (rowSize, colSize) = LA.size tab
    (constRows, lastRow) = splitAt (rowSize - 1) rows
    newTab = LA.fromRows $ constRows ++ [newRow] ++ lastRow

addGomoryCut :: Matrix R -> Vector R -> Matrix R
addGomoryCut tab gomoryRow = newTab where 
    (rowSize, colSize) = LA.size tab
    interTab = addSlackColumn tab
    rows = LA.toRows interTab
    (consRows, costRow) = splitAt (rowSize - 1) rows
    newTab = LA.fromRows (consRows ++ [gomoryRow] ++ costRow)

performGomoryCut :: ObjectiveType -> LA.Matrix R -> Int -> (R, LA.Vector R, LA.Matrix R)
performGomoryCut obj tab varIdx = (newVal, newSol, newTab) where 
    interTab = addGomoryCut tab $ getGomoryCut $ tab ! varIdx 
    newTab = updateTabDual obj interTab
    newSol = getSolution newTab
    (rowSize, colSize) = LA.size newTab
    lastVal = newTab ! (rowSize - 1) ! (colSize - 1)
    newVal = if obj == Maximization then lastVal else (-lastVal)

isClose :: R -> R -> Bool
isClose x y = diff < epsilonTol where
    diff = abs $ x - y

roundSolution :: R -> R
roundSolution num
    | diff < epsilonTol = roundCand
    | otherwise = num where
        roundCand = fromIntegral $ round num 
        diff = abs $ roundCand - num

integerSolved :: [Bool] -> [R] -> [Bool]
integerSolved = zipWith isIntSol where
    isIntSol False _ = True
    isIntSol True num = num == fromInteger (round num)

findNonIntIndex :: [Bool] -> [Bool] -> Int
findNonIntIndex intMask solMask = solIdx where
    f False solBool = True
    f True True = True
    f True False = False
    Just solIdx = elemIndex False $ zipWith f intMask solMask

getBranches :: Matrix R -> Vector R -> Int -> (Matrix R, Matrix R)
getBranches tab solVec branchIdx = (leftTab, rightTab) where
    currSolVal = solVec ! branchIdx
    colSize = LA.cols tab
    baseVec = LA.vjoin [LA.konst 0 branchIdx :: Vector R, vector [1::R], LA.konst (0::R) (colSize - branchIdx - 2)]

    leftBound = fromIntegral $ floor currSolVal
    leftRow = LA.vjoin [baseVec, vector [1::R, leftBound]]
    rightBound = fromIntegral $ ceiling currSolVal
    rightRow = LA.vjoin [-baseVec, vector [1::R, -rightBound]]

    slackTab = addSlackColumn tab
    leftTab = addNewRow slackTab leftRow
    rightTab = addNewRow slackTab rightRow

data Tree a = Nil | Node a (Tree a) (Tree a) deriving (Show)

data BranchProblem = BranchProblem {
    tableau :: Matrix R,
    solution :: Vector R,
    value :: R
} deriving (Show)


-- constructBranchAndBoundQueue :: ObjectiveType -> Vector Bool -> Vector R -> R -> [BranchProblem] -> [BranchProblem]
-- constructBranchAndBoundQueue obj intMask costVec bestVal [] = []
-- constructBranchAndBoundQueue obj intMask costVec bestVal (bp@BranchProblem{..}:bpRest)
--     | (and $ integerSolved intList currList) && (bestVal < candVal)
--     tab = tableau
--     (y, currSol, newTab) = simplexWithMixedTab obj tab
--     candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol
--     intList = LA.toList intMask 
--     currList = LA.toList currSol

-- constructBranchAndBoundBFS :: ObjectiveType -> LA.Matrix R -> LA.Vector Bool -> LA.Vector R -> R -> 

constructBranchAndBound :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Tree BranchProblem
constructBranchAndBound obj tab intMask costVec
    | infeasible = Nil
    | and $ integerSolved intList currList = currProb Nil Nil
    | otherwise = currProb leftTree rightTree where
        (y, currSol, newTab) = simplexWithMixedTab obj tab
        infeasible = y == (-infinity)
        candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol 
        currProb = Node BranchProblem {
           tableau = tab, solution = currSol, value = candVal
        }
        intList = LA.toList intMask 
        currList = LA.toList currSol
        solMask =  map (isInt . roundSolution) $ LA.toList currSol
        nextIdx = findNonIntIndex intList solMask
        (leftTab, rightTab) = getBranches tab currSol nextIdx
        leftTree = constructBranchAndBound obj leftTab intMask costVec
        rightTree = constructBranchAndBound obj rightTab intMask costVec

constructBranchAndCut :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Tree BranchProblem
constructBranchAndCut obj tab intMask costVec
    | infeasible = Nil
    | and $ integerSolved intList currList = currProb Nil Nil
    | otherwise = currProb leftTree rightTree where
        (y, currSol, newTab) = simplexWithMixedTab obj tab
        infeasible = y == (-infinity)
        candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol 
        currProb = Node BranchProblem {
           tableau = tab, solution = currSol, value = candVal
        }
        intList = LA.toList intMask 
        currList = map roundSolution $ LA.toList currSol
        solMask =  map (isInt . roundSolution) $ LA.toList currSol
        nextIdx = findNonIntIndex intList solMask
        cutTab = addGomoryCut tab $ getGomoryCut $ newTab ! nextIdx
        (leftTab, rightTab) = getBranches cutTab currSol nextIdx
        leftTree = constructBranchAndCut obj leftTab intMask costVec
        rightTree = constructBranchAndCut obj rightTab intMask costVec

constructParBranchAndCut :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Tree BranchProblem
constructParBranchAndCut obj tab intMask costVec
    | infeasible = Nil
    | and $ integerSolved intList currList = currProb Nil Nil
    | otherwise = currProb leftTree rightTree where
        (y, currSol, newTab) = simplexWithMixedTab obj tab
        infeasible = y == (-infinity)
        candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol 
        currProb = Node BranchProblem {
           tableau = tab, solution = currSol, value = candVal
        }
        intList = LA.toList intMask 
        currList = map roundSolution $ LA.toList currSol
        solMask =  map (isInt . roundSolution) $ LA.toList currSol
        nextIdx = findNonIntIndex intList solMask
        cutTab = addGomoryCut tab $ getGomoryCut $ newTab ! nextIdx
        (leftTab, rightTab) = getBranches cutTab currSol nextIdx
        (leftTree, rightTree) = runEval $ do
            leftTree <- rpar $ constructBranchAndCut obj leftTab intMask costVec
            rightTree <- rpar $ constructBranchAndCut obj rightTab intMask costVec
            _ <- rseq leftTree
            _ <- rseq rightTree
            return (leftTree, rightTree)

searchBBTreeMax :: Tree BranchProblem -> BranchProblem
searchBBTreeMax Nil = BranchProblem{tableau = LA.fromLists [[]], solution = LA.fromList [], value = -infinity}
searchBBTreeMax (Node bp Nil Nil) = bp
searchBBTreeMax (Node bp leftBpNode rightBpNode)
    | leftBpVal > rightBpVal = leftBp
    | otherwise = rightBp where
        leftBp = searchBBTreeMax leftBpNode
        rightBp = searchBBTreeMax rightBpNode
        midBpVal = value bp
        leftBpVal = value leftBp
        rightBpVal = value rightBp


branchAndBound :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> BranchProblem
branchAndBound obj tab intMask costVec = searchBBTreeMax $ constructBranchAndBound obj tab intMask costVec

branchAndCut :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> BranchProblem
branchAndCut obj tab intMask costVec = searchBBTreeMax $ constructBranchAndCut obj tab intMask costVec

obj :: ObjectiveType
obj = Maximization
intMask :: Vector Bool
intMask = LA.fromList [True, True]
costVec :: Vector R
costVec = c
xTree :: BranchProblem
xTree = searchBBTreeMax $ constructParBranchAndCut obj oldTab intMask costVec
main :: IO()
main = print xTree