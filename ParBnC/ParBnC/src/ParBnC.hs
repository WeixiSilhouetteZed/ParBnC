module ParBnC
    ( 
        ObjectiveType,
        toTableau,
        addSlackMatrix,
        constructBranchAndBound,
        constructParBranchAndBound,
        constructBranchAndCut,
        constructParBranchAndCut,
        searchBBTreeMax
    ) where

import Numeric.IEEE ( IEEE(infinity) )
import Data.List ( elemIndex )
import Numeric.LinearAlgebra as LA
    ( (<.>),
      dropColumns,
      fromLists,
      takeColumns,
      cols,
      flatten,
      fromColumns,
      fromRows,
      rows,
      toColumns,
      toRows,
      cmap,
      find,
      maxElement,
      maxIndex,
      scalar,
      sumElements,
      diagl,
      size,
      vector,
      subVector,
      takesV,
      toList,
      vjoin,
      fromList,
      Matrix,
      Konst(konst),
      Indexable((!)),
      R,
      Vector )
import Control.Parallel(par, pseq)
import Control.Parallel.Strategies
    ( rdeepseq, rparWith, runEval, NFData )
import Control.Monad ( when )
import Control.DeepSeq ( NFData(..) )

data MatConstraints = MatVec [[Double]] [Double] deriving Show

data VariableType = INTEGER | CONTINUOUS deriving (Show, Eq)

data ProblemType = LP | MIP deriving Show

data ObjectiveType = Maximization | Minimization deriving (Show, Eq)

epsilonTol :: R
epsilonTol = 1e-10

failNodeThreshold :: Int
failNodeThreshold = 100

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
isImprovable obj tab = any (costCheck obj . roundSolution) $ LA.toList cost where 
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
        | elem == 0 = infinity::R
        | val < 0 = infinity::R
        | otherwise = val
        where 
            elem = rowEq ! column
            lastColVal = rowEq ! (LA.size rowEq - 1)
            val =  lastColVal / (rowEq ! column)
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
    rowSize = LA.rows tab
    newTableau = LA.fromRows $ zipWith updateRow [0..(rowSize - 1)] $ LA.toRows tab

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

updateTab :: ObjectiveType -> LA.Matrix R -> Int -> LA.Matrix R
updateTab obj tab counter
    | not $ isImprovable obj tab = tab
    | counter == 0 = LA.fromLists [[]]
    | otherwise = updateTab obj newTab (counter - 1) where 
        pivotPos = getPivotPosition obj tab
        newTab = pivotStep tab pivotPos


updateTabDebug :: ObjectiveType -> LA.Matrix R -> Int -> LA.Matrix R
updateTabDebug obj tab counter
    | counter == 0 = tab
    | not $ isImprovable obj tab = tab
    | otherwise = updateTabDebug obj newTab (counter - 1) where 
        pivotPos = getPivotPosition obj tab 
        newTab = pivotStep tab pivotPos

updateTabDual :: ObjectiveType -> LA.Matrix R -> LA.Matrix R
updateTabDual obj tab 
    | not $ isImprovableDual obj tab = tab
    | otherwise = updateTabDual obj newTab where 
        pivotPos = getPivotPositionDual obj tab
        newTab = pivotStep tab pivotPos

isImprovableMixed :: LA.Matrix R -> Bool
isImprovableMixed tab = any ((> 0) . roundSolution) (LA.toList costs) where
    lastRow = last $ LA.toRows tab
    costs = subVector 0 (LA.cols tab - 1) lastRow

getPhaseOneTab :: LA.Matrix R -> (LA.Matrix R, LA.Vector R, Bool)
getPhaseOneTab tab = (newTab, oldCost, needed) where 
    (rowSize, colSize) = LA.size tab
    (constRows, [oldCost]) = splitAt (rowSize - 1) $ LA.toRows tab
    updateRow :: LA.Vector R -> LA.Vector R
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
            lastColVal = rowEq ! (LA.size rowEq - 1)
            val =  lastColVal / (rowEq ! column)
    restrictions = map getElem $ init (LA.toRows tab)
    minRatio = minimum restrictions
    row =  last [idx | (idx, ratio) <- zip [0..colSize] restrictions, ratio == minRatio]

updateMixedTab :: LA.Matrix R -> Int -> LA.Matrix R
updateMixedTab tab counter
    | not $ isImprovableMixed tab = tab
    | counter == 0 = LA.fromLists [[]]
    | otherwise = updateMixedTab newTab (counter - 1) where 
        pivotPos = getPivotPositionMixed tab
        newTab = pivotStep tab pivotPos

simplexWithTab :: ObjectiveType -> LA.Matrix R -> (R, LA.Vector R, LA.Matrix R)
simplexWithTab obj tab = (optVal, solution, lastTab) where 
    lastTab = updateTab obj tab failNodeThreshold
    solution = getSolution lastTab
    (rowSize, colSize) = LA.size lastTab
    lastVal = lastTab ! (rowSize - 1) ! (colSize - 1)
    optVal = if obj == Maximization then lastVal else (-lastVal)

simplexWithMixedTab :: ObjectiveType -> LA.Matrix R -> (R, LA.Vector R, LA.Matrix R)
simplexWithMixedTab obj tab
    | failedInter = (-infinity, getSolution phaseOneTab, phaseOneTab)
    | infeasible = (-infinity, getSolution interTab, interTab)
    | failedNode = (-infinity, getSolution interTab, interTab)
    | otherwise  = (lastVal, solution, lastTab) where
        (phaseOneTab, oldCost, needed) = getPhaseOneTab tab 
        interTab
            | needed = updateMixedTab phaseOneTab failNodeThreshold
            | otherwise = phaseOneTab
        failedInter = LA.fromLists [[]] == interTab
        (rowSize, colSize) = LA.size interTab
        lastPhaseOneVal = interTab ! (rowSize - 1) ! (colSize - 1)
        infeasible = not $ isClose lastPhaseOneVal (0::R)

        (constRows, _) = splitAt (rowSize - 1) $ LA.toRows interTab

        phaseTwoTab = LA.fromRows $ constRows ++ [oldCost]
        lastTab = updateTab obj phaseTwoTab failNodeThreshold
        failedNode = lastTab == LA.fromLists [[]]
        solution = getSolution lastTab
        lastVal = lastTab ! (rowSize - 1) ! (colSize - 1)
        optVal = lastVal

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

addSlackMatrix :: LA.Matrix R -> LA.Matrix R
addSlackMatrix tab = newTab where
    (rowSize, colSize) = LA.size tab
    slackIdentity = LA.diagl $ replicate (rowSize - 1) (1::R)
    zeroVec = LA.konst (0::R)  (rowSize - 1)
    newSlackMat = LA.fromRows $ LA.toRows slackIdentity ++ [zeroVec]
    (prevCols, lastCol) = splitAt (colSize - 1) $ LA.toColumns tab
    newTab = LA.fromColumns (prevCols ++ LA.toColumns newSlackMat ++ lastCol)

isClose :: R -> R -> Bool
isClose x y = diff < epsilonTol where
    diff = abs $ x - y

roundSolution :: R -> R
roundSolution num
    | diff < epsilonTol = roundCand
    | otherwise = num where
        roundCand = fromIntegral $ round num 
        diff = abs $ roundCand - num

integerSolved :: [Bool] -> [Bool] -> [Bool]
integerSolved = zipWith isIntSol where
    isIntSol False _ = True
    isIntSol True mask = mask

fromJust :: Maybe a -> a
fromJust (Just a) = a
fromJust Nothing = error "non-existing index"

findNonIntIndex :: [Bool] -> [Bool] -> Int
findNonIntIndex intMask solMask = solIdx where
    f False solBool = True
    f True True = True
    f True False = False
    maybeSolIdx = elemIndex False $ zipWith f intMask solMask
    solIdx = fromJust maybeSolIdx
        

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
instance NFData a => NFData (Tree a) where
    rnf Nil = ()
    rnf (Node l a r) = rnf l `seq` rnf a `seq` rnf r

data BranchProblem = BranchProblem {
    solution:: Vector R, value :: R
} deriving (Show)
instance NFData BranchProblem where
    rnf BranchProblem {solution = s, value = v} = rnf s `seq` rnf v

constructBranchAndBound :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Int -> Tree BranchProblem
constructBranchAndBound obj tab intMask costVec maxDepth
    | infeasible = Nil
    | maxDepth == 0 = Nil
    | and $ integerSolved intList solMask = currProb Nil Nil
    | otherwise = currProb leftTree rightTree where
        (y, currSol, newTab) = simplexWithMixedTab obj tab
        infeasible = y == (-infinity)
        candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol 
        currProb = Node BranchProblem {
           solution = currSol, value = candVal
        }
        intList = LA.toList intMask 
        solMask =  map (isInt . roundSolution) $ LA.toList currSol
        nextIdx = findNonIntIndex intList solMask
        (leftTab, rightTab) = getBranches tab currSol nextIdx
        leftTree = constructBranchAndBound obj leftTab intMask costVec (maxDepth - 1)
        rightTree = constructBranchAndBound obj rightTab intMask costVec (maxDepth - 1)

constructParBranchAndBound :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Int -> Int -> Tree BranchProblem
constructParBranchAndBound obj tab intMask costVec maxDepth parDepth
    | infeasible = Nil
    | maxDepth == 0 = Nil
    | and $ integerSolved intList solMask = currProb Nil Nil
    | otherwise = currProb leftTree rightTree where
        (y, currSol, newTab) = simplexWithMixedTab obj tab
        infeasible = y == (-infinity)
        candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol 
        currProb = Node BranchProblem {
            solution = currSol, value = candVal
        }
        intList = LA.toList intMask 
        solMask =  map (isInt . roundSolution) $ LA.toList currSol
        nextIdx = findNonIntIndex intList solMask
        (leftTab, rightTab) = getBranches tab currSol nextIdx
        (leftTree, rightTree) 
            | parDepth == 0 = (leftRawTree, rightRawTree)
            | otherwise = leftTree `par` rightTree `pseq` (leftTree, rightTree) where 
                leftRawTree = constructBranchAndBound obj leftTab intMask costVec (maxDepth - 1)
                rightRawTree = constructBranchAndBound obj rightTab intMask costVec (maxDepth - 1)
                leftTree = constructParBranchAndBound obj leftTab intMask costVec (maxDepth - 1) (parDepth - 1)
                rightTree = constructParBranchAndBound obj rightTab intMask costVec (maxDepth - 1) (parDepth - 1)

constructBranchAndCut :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Int -> Tree BranchProblem
constructBranchAndCut obj tab intMask costVec maxDepth
    | maxDepth == 0 = Nil
    | infeasible = Nil
    | and $ integerSolved intList solMask = currProb Nil Nil
    | otherwise = currProb leftTree rightTree where
        (y, currSol, newTab) = simplexWithMixedTab obj tab
        infeasible = y == (-infinity)
        candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol 
        currProb = Node BranchProblem {
           solution = currSol, value = candVal
        }
        intList = LA.toList intMask 
        solMask =  map (isInt . roundSolution) $ LA.toList currSol
        nextIdx = findNonIntIndex intList solMask
        cutTab = addGomoryCut tab $ getGomoryCut $ newTab ! nextIdx
        (leftTab, rightTab) = getBranches cutTab currSol nextIdx
        leftTree = constructBranchAndCut obj leftTab intMask costVec (maxDepth - 1)
        rightTree = constructBranchAndCut obj rightTab intMask costVec (maxDepth - 1)

constructParBranchAndCut :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Int -> Int -> Tree BranchProblem
constructParBranchAndCut obj tab intMask costVec maxDepth parDepth
    | infeasible = Nil
    | maxDepth == 0 = Nil
    | and $ integerSolved intList solMask = currProb Nil Nil
    | otherwise = currProb leftTree rightTree where
        (y, currSol, newTab) = simplexWithMixedTab obj tab
        infeasible = y == (-infinity)
        candVal = costVec <.> LA.subVector 0 (LA.size costVec) currSol 
        currProb = Node BranchProblem {
            solution = currSol, value = candVal
        }
        intList = LA.toList intMask 
        solMask =  map (isInt . roundSolution) $ LA.toList currSol
        nextIdx = findNonIntIndex intList solMask
        cutTab = addGomoryCut tab $ getGomoryCut $ newTab ! nextIdx
        (leftTab, rightTab) = getBranches cutTab currSol nextIdx
        (leftTree, rightTree) = runEval $ do
            if parDepth == 0 then do
                leftRawTree <- rparWith rdeepseq $ constructBranchAndCut obj leftTab intMask costVec (maxDepth - 1)
                rightRawTree <- rparWith rdeepseq $  constructBranchAndCut  obj rightTab intMask costVec (maxDepth - 1)
                return (leftRawTree, rightRawTree)
            else do 
                leftTree <- rparWith rdeepseq  $  constructParBranchAndCut  obj leftTab intMask costVec (maxDepth - 1) (parDepth - 1)
                rightTree <- rparWith rdeepseq  $  constructParBranchAndCut  obj rightTab intMask costVec (maxDepth - 1) (parDepth - 1)
                return (leftTree, rightTree)

searchBBTreeMax :: Tree BranchProblem -> BranchProblem
searchBBTreeMax Nil = BranchProblem{solution = LA.fromList [], value = -infinity}
searchBBTreeMax (Node bp Nil Nil) = bp
searchBBTreeMax (Node bp leftBpNode Nil) = searchBBTreeMax leftBpNode
searchBBTreeMax (Node bp Nil rightBpNode) = searchBBTreeMax rightBpNode
searchBBTreeMax (Node bp leftBpNode rightBpNode)
    | leftBpVal > rightBpVal = leftBp
    | otherwise = rightBp where
        leftBp = searchBBTreeMax leftBpNode
        rightBp = searchBBTreeMax rightBpNode
        midBpVal = value bp
        leftBpVal = value leftBp
        rightBpVal = value rightBp

