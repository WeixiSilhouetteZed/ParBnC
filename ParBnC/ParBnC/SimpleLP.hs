
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

constructParBranchAndBound :: ObjectiveType -> Matrix R -> Vector Bool -> Vector R -> Tree BranchProblem
constructParBranchAndBound obj tab intMask costVec
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
        (leftTree, rightTree) = runEval $ do
            leftTree <- rpar $ constructParBranchAndBound obj leftTab intMask costVec
            rightTree <- rpar $ constructParBranchAndBound obj rightTab intMask costVec
            _ <- rseq leftTree
            _ <- rseq rightTree
            return (leftTree, rightTree)

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
            leftTree <- rpar $ constructParBranchAndCut obj leftTab intMask costVec
            rightTree <- rpar $ constructParBranchAndCut obj rightTab intMask costVec
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

testMat :: LA.Matrix R
testMat = LA.fromLists [[ 8.55802933e+00,  5.08098357e+00,  4.05670235e+00,
          1.70539448e-01,  9.88113522e-01,  9.51424163e+00,
          4.65492716e+00,  1.90727612e+00,  5.37358570e+00,
          1.36738550e+00,  6.40485564e+00,  3.71785512e+00,
          8.13439689e+00,  2.68577818e+00,  9.18271510e-01,
          4.30410666e+00,  4.51763462e+00,  2.44787851e+00,
         -1.05393211e+00,  5.02100056e+00,  2.70590292e+00,
          4.89202969e+00,  2.56371258e-01,  1.04199959e+00,
          9.13579562e-01,  4.33041124e+00, -5.73797420e-01,
          1.43787363e+01, -4.38927210e-01,  1.07574075e+01,
          6.94420200e+00,  1.21148648e+00,  2.05357889e+00,
          3.91646596e+00,  7.71444922e+00,  4.81460905e+00,
          4.97915607e+00, -1.32186368e+00,  7.17339004e+00,
          0.00000000e+00,  1.24477320e+00],
        [ 3.78973252e+00,  5.39165567e+00,  4.73560311e+00,
         -2.99648058e+00,  1.21962935e+00,  8.30267839e+00,
         -1.69175460e+00,  5.94715567e+00,  3.52048334e+00,
         -1.83081646e+00,  8.57222917e-02,  2.52033707e+00,
          7.93763360e-01, -2.18124212e+00,  3.58719636e-01,
          1.87471386e+00,  3.71926286e-01,  1.28998283e+00,
          5.34280136e-01,  1.62150919e+00, -1.14870070e-01,
          6.07361496e+00, -6.02154260e-01, -7.69697740e-01,
         -4.39963840e+00,  5.17414125e+00, -2.29737403e+00,
          1.09790220e+01,  2.74070256e+00,  6.76321042e+00,
          7.35487895e+00, -7.03198150e-01,  5.73927977e+00,
          2.55686961e+00, -1.35041980e-01,  4.40900597e+00,
          4.82178956e+00, -2.23491972e+00,  4.86067783e+00,
          2.40425748e+00,  1.58275704e+00],
        [ 0.00000000e+00,  1.41092985e+00, -5.83174639e+00,
          0.00000000e+00,  1.89873510e+00,  7.35096780e+00,
          0.00000000e+00,  4.96678709e+00,  0.00000000e+00,
          1.12067527e+00,  2.96025607e+00,  2.92596715e+00,
          6.30798676e-01,  0.00000000e+00,  0.00000000e+00,
         -1.21091552e+00,  1.34618664e+00, -1.06612044e+00,
          0.00000000e+00,  3.40515967e+00,  3.90245858e-01,
          0.00000000e+00,  2.11528773e+00,  0.00000000e+00,
          1.85856637e+00,  4.68158213e+00,  0.00000000e+00,
          0.00000000e+00,  4.71983881e+00,  8.09808747e+00,
          1.59856030e+00,  7.84317986e-01,  4.75098200e+00,
          0.00000000e+00,  4.43264372e+00,  5.41917616e+00,
          0.00000000e+00,  4.87942234e-01,  0.00000000e+00,
          0.00000000e+00,  1.65177997e+00],
        [ 2.61009732e-02, -1.53319273e+00,  3.76155269e+00,
          1.39991615e-01,  3.73240039e+00,  9.36258830e+00,
         -2.88935804e+00,  3.34919519e+00,  4.44405953e+00,
          7.74013269e-01,  1.53816014e-01,  2.89742203e+00,
          4.86253622e+00,  3.51039710e+00,  2.47898082e+00,
         -1.45894670e-01,  1.84304770e+00,  1.43060423e+00,
         -3.08094998e+00,  5.10193660e+00,  2.41547721e+00,
          4.55272871e+00, -4.02905560e-01,  1.12917680e+00,
         -1.48076205e+00,  5.64660532e+00,  4.04365041e+00,
          5.93020804e+00,  2.27512066e+00,  1.46313810e+01,
          6.67304000e+00,  1.68279555e+00,  4.46913224e+00,
          6.20236297e+00, -1.56606150e-02,  3.60558090e+00,
          3.07739087e-01, -7.70373260e-01, -2.15885965e+00,
          1.34921718e+00,  1.66451727e+00],
        [ 6.01898100e+00,  6.38419215e+00, -3.31982976e+00,
         -3.67756710e+00, -2.19743009e+00,  6.88050384e+00,
         -7.57573520e-01,  4.00546692e+00,  3.01619234e+00,
         -2.33840150e-01,  2.00533125e-01,  1.87835523e+00,
          2.82433526e+00,  1.36647525e+00,  0.00000000e+00,
          1.19259788e+00,  3.99640537e+00,  2.40557870e+00,
          6.92127408e-01,  3.02196894e+00,  1.83266552e+00,
          1.69726175e+00,  5.80541990e-01, -2.49704960e-01,
          2.69149047e+00,  3.73953825e+00, -8.65052180e-01,
          0.00000000e+00,  1.11821293e+00,  1.87636991e+00,
          2.24524046e+00,  2.90274209e+00,  2.15526954e+00,
          7.10948923e-01, -3.79028650e-01,  4.36108158e+00,
          7.62014220e-01,  7.73967388e-01, -5.21616310e-01,
          1.43582727e+00, -1.50905342e+00],
        [ 1.67592609e+00,  3.90734139e+00,  3.93340908e+00,
          1.84987556e+00,  3.66194615e+00,  8.39590385e+00,
          1.74315589e+00,  2.89445263e+00,  2.96956282e+00,
          1.39904545e+00,  5.39544089e+00,  3.71394558e+00,
          5.59628241e+00,  4.11578695e-01,  0.00000000e+00,
          2.35379463e+00,  4.46160665e+00,  3.11378221e+00,
          2.60688675e+00,  4.11131053e+00,  1.88239241e-01,
          4.79412800e+00, -9.63843720e-01,  6.60047806e-01,
          1.45076345e+00,  1.35828216e+00,  5.63101510e+00,
          1.79175313e+01,  5.50366416e-01,  1.38703131e+01,
          2.85710494e+00,  1.14559899e+00,  1.32767901e+00,
          4.69238799e+00, -2.89656580e-02,  1.02390982e+00,
          6.12918812e+00,  6.32498711e+00,  2.61522296e+00,
          7.38642325e+00, -2.40656148e+00],
        [ 1.99440062e+00,  2.24594854e+00,  3.15265925e+00,
         -2.69815100e-01,  3.09091650e+00,  8.02751680e+00,
          2.03668531e+00,  4.47673203e+00,  1.39165703e-01,
          1.83140285e+00,  1.58460758e+00,  3.80736718e+00,
          4.89144477e+00,  4.11759502e+00,  3.58672394e+00,
          3.22287026e-01,  5.98961579e+00,  4.70625936e+00,
          3.18991164e+00,  2.13751483e+00,  2.55254222e+00,
          2.38313484e+00,  3.94302033e+00, -1.99982820e-01,
          4.62241025e+00,  4.31615871e+00, -1.91048203e+00,
          6.62304625e+00,  1.00157403e+00,  5.90643376e+00,
          5.33580084e+00, -1.28622059e+00,  4.39283976e+00,
         -6.64619450e-01, -1.19535120e-01, -4.50226920e-01,
          3.68035852e+00, -3.05749993e+00,  2.89365831e+00,
          0.00000000e+00,  5.78079509e+00],
        [ 1.70541190e-01,  1.05672605e+00,  9.32983159e-01,
         -2.51816886e+00,  2.56527588e+00,  2.39529598e+00,
         -2.69249025e+00,  9.46569169e-01,  5.28949508e+00,
          2.77355103e+00,  1.69767564e+00, -6.60548010e-01,
          2.81550432e+00,  9.66242985e-01,  0.00000000e+00,
          2.22605161e+00, -1.52767135e+00,  1.30803726e-01,
         -9.30562870e-01, -1.26493276e+00,  1.39442071e+00,
          1.86190713e+00, -8.51187200e-01, -2.36222267e+00,
          3.68023485e+00,  8.09864922e+00,  2.95775326e+00,
          0.00000000e+00,  4.77913648e-01,  3.84891171e+00,
          5.82619616e+00,  4.49969817e-01,  4.98677206e+00,
         -1.68430195e+00, -1.90074500e+00,  3.73359398e+00,
          9.89424777e-01,  1.59632148e+00,  3.56740158e+00,
          4.57276879e+00,  4.82634712e+00],
        [-1.37415977e+00,  3.57615240e+00, -2.98127568e+00,
          2.25709639e+00, -2.42255410e-02,  8.46147921e+00,
          3.13363926e+00,  3.62309601e+00,  3.71438829e+00,
          3.47784378e+00,  6.02265146e-01,  5.26335461e+00,
          5.04808069e+00,  1.56791037e+00,  3.49816346e+00,
          2.73396166e+00,  3.63194642e+00,  3.92068246e+00,
          1.48898244e+00,  4.57467495e+00,  2.00659112e-01,
          4.74400200e+00,  3.67949722e+00,  2.67850121e+00,
          2.05369331e+00,  6.03105826e+00,  4.58958971e+00,
          0.00000000e+00, -8.19212690e-02,  9.28263548e+00,
          1.94527753e+00, -1.41046178e+00,  2.25202180e+00,
          1.57773361e+00, -2.35546262e+00,  6.99870359e-01,
          3.17875761e+00, -2.35109230e-01,  1.76944908e+00,
          1.44342648e+00,  3.17735267e+00],
        [-3.89155930e-02, -2.86479770e-01,  1.04674493e+00,
          2.84609931e+00, -8.92290930e-01,  6.62712196e+00,
         -7.47669171e+00,  5.43758111e-01,  4.09697017e-01,
         -6.17448210e-01,  2.84221716e+00,  4.66251741e+00,
          2.52663166e+00,  2.56339014e-01, -4.51192050e-01,
          1.46925385e+00,  3.26932922e+00, -1.43059838e+00,
          1.54181660e+00,  2.10092591e-01, -4.56891639e+00,
         -4.03314918e+00,  5.33411913e-01, -6.28345526e+00,
         -1.67081579e+00,  3.17127817e+00,  2.65700626e+00,
          1.37790984e+01,  4.82744108e+00,  5.10137271e+00,
          5.23418632e+00,  8.80021171e-01,  2.93259093e+00,
         -4.22439040e-01, -2.40646241e+00,  2.66149864e+00,
          5.49985814e+00,  1.64507503e+00,  3.67902059e+00,
          2.55129937e+00,  2.69700102e+00],
        [ 1.13272479e+00,  2.33524930e+00,  1.72358699e+00,
          1.44537388e+00, -1.06168456e+00,  8.15596477e+00,
          3.38970412e+00,  2.11380009e+00,  9.81903098e-01,
         -2.83433045e+00,  7.16445864e-02, -2.50454439e+00,
          4.08541408e+00, -8.55609580e-01,  1.02439976e+00,
         -8.81429170e-01,  6.10201633e+00,  3.21134970e-03,
          2.02178436e+00,  2.13179228e+00,  5.75846408e-01,
          2.27214705e+00,  6.56244912e-01, -3.53354010e-01,
          1.99693380e+00,  8.23127935e+00,  4.31138733e+00,
          8.86213005e+00,  8.40202693e-01,  1.09949976e+01,
          5.24879274e+00,  1.58525049e+00,  2.45242411e+00,
          3.07694599e-01,  3.93565132e+00,  2.25083252e+00,
          3.89677920e+00,  3.76765909e+00, -4.44259090e-01,
          4.73642005e+00,  2.38688417e-01],
        [ 0.00000000e+00, -4.25605140e-01,  2.06540945e+00,
          0.00000000e+00,  1.55603751e+00,  0.00000000e+00,
          4.08015456e+00,  3.27094854e+00,  0.00000000e+00,
          0.00000000e+00,  3.69530781e+00, -1.42422911e+00,
          5.63769058e+00,  1.40163932e-01,  0.00000000e+00,
          3.22098223e-01,  4.96335519e+00,  1.38614827e+00,
         -1.14042200e+00,  1.26291993e+00,  5.39424804e-01,
          4.36448648e+00, -2.94676335e+00, -6.20702670e-01,
          1.31907802e+00,  2.76059478e+00,  0.00000000e+00,
          7.58887753e+00,  0.00000000e+00,  0.00000000e+00,
          0.00000000e+00,  0.00000000e+00,  5.38331509e+00,
         -1.08476863e+00, -1.10476950e-01, -1.07329097e+00,
         -6.83200100e-01, -3.39417439e+00,  1.15939498e+00,
          1.96903748e+00,  0.00000000e+00],
        [ 3.00148428e+00,  3.03791772e+00,  1.13044651e+00,
          5.98270167e-01, -5.88513970e-01,  5.58219459e+00,
         -5.18460120e-01,  6.89159197e+00,  5.75212099e+00,
          1.54374436e+00,  1.97122332e+00,  9.47314541e-01,
          4.93600366e+00, -3.86845030e-01,  0.00000000e+00,
          2.22885137e+00,  4.11049191e+00,  2.54173626e+00,
         -5.80616580e-01, -1.41245533e+00,  1.02608933e+00,
          2.12035325e-01,  4.08187892e+00, -4.85527059e+00,
          1.07867437e+00,  2.60518461e+00,  3.53469207e-01,
          0.00000000e+00,  1.79036683e-01,  7.57875590e+00,
          3.17054183e+00,  1.06735825e-01,  4.29574550e+00,
         -2.22637857e+00,  1.18491129e+00, -1.14106160e-01,
          6.03291015e+00, -3.00830481e+00, -1.49327178e+00,
         -3.75853295e+00,  2.20375682e-01],
        [ 2.94807524e+00,  7.36483017e+00,  3.30342900e+00,
         -5.00278340e+00,  4.02847028e+00,  6.81477225e+00,
         -3.59855430e-01, -1.04279480e+00,  1.85879354e+00,
          3.07612404e+00,  1.74689016e+00,  3.24554513e+00,
          2.84158747e+00, -1.46776806e+00,  0.00000000e+00,
         -1.81208880e+00,  6.44058926e+00,  8.11696060e-01,
         -2.82907126e+00,  5.74489381e+00,  0.00000000e+00,
          2.98773410e+00, -1.52032845e+00,  8.70336946e-01,
         -1.31448829e+00,  4.84003525e+00,  3.09545099e+00,
          8.25773821e+00,  3.94985144e+00,  0.00000000e+00,
          2.67635643e+00,  1.71790000e+00,  2.29919071e+00,
         -2.46758827e+00,  3.99821892e+00,  2.00351914e+00,
          2.18080034e+00,  2.71995396e+00,  2.34302001e+00,
         -2.95378905e+00,  2.87717327e+00],
        [ 2.94352139e+00,  6.91425653e+00, -2.99080216e+00,
          1.20478939e+00,  1.90049733e+00,  9.44779712e+00,
          3.27138881e+00,  5.21444430e+00,  4.10987261e+00,
          4.73832757e+00,  3.80294147e+00,  6.61923928e-01,
          4.89463966e+00, -1.28452088e+00,  2.10457394e+00,
          2.95025976e+00,  1.15073922e+00,  2.43254289e+00,
          2.20495754e+00,  9.65852037e-01, -1.03450334e+00,
          7.35583243e+00,  4.83275904e+00, -5.21309244e+00,
          2.39772120e+00,  2.88808154e+00,  3.87673462e+00,
          7.31467417e+00,  3.18296908e+00,  8.56942054e+00,
          7.61133048e+00,  1.70550607e+00,  2.95255244e+00,
         -1.12203484e+00, -8.66945440e-02,  2.35081814e+00,
          1.84806136e+00, -2.78541870e-01,  3.38853900e+00,
          2.25566255e+00,  6.22611510e+00],
        [ 6.58092482e+00,  5.84640955e+00,  6.56666307e+00,
          6.86444343e+00,  2.36867155e+00,  1.23597871e+01,
          6.27614950e+00,  3.58121971e+00,  3.86671186e+00,
         -1.44626480e+00,  4.10278775e+00,  3.74692051e+00,
          5.62358060e+00,  1.99712281e+00,  6.26862335e+00,
          4.84295163e+00,  6.07887022e+00,  4.48699607e+00,
          1.76792203e+00,  2.60763765e+00,  1.18295172e+01,
          3.37341242e+00,  3.16272017e+00,  2.98139752e+00,
          5.54549465e+00,  7.36785634e+00,  1.30435983e+00,
          1.74445311e+01,  8.06935409e+00,  1.10053901e+01,
          9.85310597e+00,  1.38597839e+00,  6.58295927e+00,
         -5.01773510e-01,  2.84158077e+00,  7.29752225e+00,
          6.28441927e+00,  2.77734352e+00,  2.24841079e+00,
          0.00000000e+00,  3.28039107e+00],
        [ 2.42098187e+00,  5.79505507e+00, -2.44200470e-01,
          3.55766427e+00,  2.69428559e+00,  1.09950103e+01,
          4.92907723e+00,  1.20376676e+01,  8.16093090e+00,
          1.63558128e+00,  2.28334022e+00,  3.15555343e+00,
          1.43763678e+00, -1.11882925e+00,  0.00000000e+00,
          2.76708918e+00,  1.09354213e+00,  6.17562172e+00,
          5.98389361e+00,  1.06366983e+01,  1.99272939e+00,
          5.14680400e+00,  5.07554913e-01,  3.22439819e+00,
          4.06792525e+00,  8.28662859e+00,  2.05039190e+00,
          1.51808510e+01,  3.39908052e+00,  1.07511005e+01,
          6.92577471e+00,  4.93629450e+00,  5.66066169e+00,
          3.44410041e+00,  4.22277850e+00,  2.82936615e+00,
          2.89310073e+00,  1.27490481e+00,  5.04159107e+00,
          7.02869522e+00,  4.03603737e+00],
        [ 1.33278520e+00,  2.88285438e+00,  1.13517170e+00,
         -1.24314995e+00,  1.87342058e+00,  6.33845264e+00,
          2.34621193e+00,  6.01207815e+00,  5.83646637e+00,
         -9.13607830e-01,  3.89309436e+00,  1.42053807e+00,
          5.95836318e+00, -4.20497760e-01,  2.79107139e+00,
          4.94876463e-01,  4.99016829e+00,  2.65419964e+00,
          4.17465223e+00,  1.49464731e+00,  3.98248600e+00,
          2.79359369e+00,  6.40513771e-01, -8.88397180e-01,
          9.78265063e-01,  4.25443838e+00, -1.24791683e+00,
          9.54365600e+00,  2.65899589e+00,  6.79009952e+00,
          1.71624815e+00, -2.27659482e+00,  4.70902001e+00,
          1.04744205e+00, -2.42286400e-01,  7.12190169e-01,
          3.15634749e+00,  1.07800935e+00,  4.88937227e-01,
          1.50897419e+00,  6.88403550e-01],
        [ 1.63677261e+00,  4.79746106e+00,  5.27942055e+00,
          8.23781644e-01,  6.09362445e+00,  1.01323883e+01,
          1.83229707e+00,  5.21626802e+00,  4.62966883e+00,
          3.81801579e-01,  4.14689374e+00,  5.43608504e-01,
          1.96137770e+00,  6.81642998e-01,  0.00000000e+00,
          2.90237894e+00,  4.23014077e+00,  1.73404543e+00,
          4.31871466e+00,  4.86342895e-01,  4.20864641e+00,
          1.76147609e+00,  6.97133573e-01, -4.99077060e-01,
          1.95055516e+00,  1.34855300e+00, -8.67194170e-01,
          9.29767133e+00,  1.10319112e+00,  7.57949711e+00,
         -7.60037140e-02, -2.41747851e+00,  7.62045476e+00,
          4.57073563e-01,  1.19772756e+00,  2.71126268e+00,
          5.74171103e+00, -1.09330615e+00,  2.37814460e+00,
          1.28970118e-02,  8.10887930e-01],
        [ 6.26159818e+00, -1.35305883e+00,  2.86954936e+00,
         -1.11316892e+00,  2.57661713e+00,  9.30933498e+00,
          3.06308544e+00,  5.08427639e+00,  1.89127092e+00,
         -1.70452159e+00,  4.20115577e+00,  5.79003093e+00,
          7.47642237e+00, -1.64436009e+00,  5.16045257e-01,
          9.39873546e-01,  2.77632762e+00, -2.12092840e-01,
          2.58148788e+00,  1.25339073e+00,  1.32689080e+00,
          3.80874292e+00,  2.06755565e+00, -2.91063635e+00,
         -1.93720262e+00,  5.49350070e+00, -1.68757324e+00,
          5.97808241e+00,  4.82110141e+00,  9.79326788e+00,
          6.10731466e+00, -1.83295487e+00,  2.31634907e+00,
         -1.68043778e+00, -4.25102350e-01,  2.60087560e-01,
          4.62675999e-01,  5.54072369e+00, -2.50497145e+00,
          0.00000000e+00,  7.07081584e+00],
        [ 1.53103807e+00,  6.68285622e+00,  1.20915616e+00,
          2.18860238e-01,  1.84858458e+00,  8.74018739e+00,
          1.26774001e+00,  3.80042195e+00,  4.33482229e+00,
          1.69137113e+00,  5.31847126e+00,  1.93020110e+00,
          2.24045161e+00, -2.04877704e+00,  0.00000000e+00,
          1.64375443e+00,  3.00937365e+00,  1.86003920e-01,
          3.22600936e+00,  5.50517076e+00,  1.93628993e+00,
         -1.47534308e+00,  2.76773930e+00, -2.49555343e+00,
         -1.19115859e+00,  2.61393756e+00,  1.81546703e+00,
          1.02256684e+01,  1.54544263e-01,  6.45886964e+00,
          4.16326535e+00, -7.10132360e-01,  1.39965561e+00,
          1.22539214e-02,  1.51200990e+00,  7.19721762e-01,
          4.72972771e+00,  6.34568579e-01,  2.86324269e+00,
          0.00000000e+00,  2.71691708e+00],
        [ 1.76129815e-01,  6.36072579e+00, -2.57552376e+00,
          2.80830707e+00,  3.12544342e+00,  1.09093216e+01,
         -2.87349620e-01,  1.73167646e+00,  2.98249396e+00,
          2.73160740e+00,  2.69935711e+00,  4.90879739e+00,
          3.91893326e+00,  2.50774015e+00,  3.62919458e+00,
          7.01820576e-01,  1.62188437e+00,  2.48955909e+00,
         -2.19955154e+00,  3.02913442e+00,  2.28357810e+00,
          3.07257535e+00,  3.50391185e+00, -9.44399120e-01,
          1.79885839e+00,  3.12091534e+00, -4.28101787e+00,
          1.24932327e+01, -2.17897474e+00,  6.22784587e+00,
          2.45275604e+00,  6.61851628e-01,  7.26574837e-02,
         -2.52318930e+00, -1.21978588e+00,  2.65185319e-01,
          4.70462015e+00,  2.05305571e+00, -3.69469790e+00,
         -3.96002400e-01, -3.74756940e-01],
        [ 2.70341105e+00,  5.82605308e+00,  4.20057169e+00,
          0.00000000e+00, -3.54330874e+00,  7.72547694e+00,
          1.69547060e+00, -1.89262432e+00,  0.00000000e+00,
          0.00000000e+00, -1.27671958e+00,  0.00000000e+00,
          8.02305566e+00, -8.50978200e-01,  0.00000000e+00,
         -9.74770010e-01,  2.00475423e+00,  0.00000000e+00,
          0.00000000e+00,  2.72024250e-01,  0.00000000e+00,
          1.21011253e-01, -1.95060353e+00,  0.00000000e+00,
          1.06741777e+00,  7.70752980e+00,  0.00000000e+00,
          1.20503782e+01,  6.53804217e+00,  8.01607270e+00,
          4.84606852e+00,  0.00000000e+00,  4.69340893e+00,
          2.19179302e+00,  0.00000000e+00,  0.00000000e+00,
          4.50799510e+00, -2.03897369e+00,  0.00000000e+00,
          2.36756390e+00, -1.48917759e+00],
        [ 2.48796747e+00,  3.13671913e+00,  6.02843149e-02,
          1.90283977e+00,  1.27426939e+00,  1.00143627e+01,
          1.73583399e+00,  1.40487244e+00,  1.46205593e+00,
          5.26250089e+00,  3.65953797e+00,  7.60058547e+00,
          3.45415471e+00,  1.15406535e+00, -1.54843937e+00,
          4.16676136e+00, -1.88535900e-02, -1.00512727e+00,
          0.00000000e+00,  3.21203010e+00,  3.50782466e-01,
          1.18845490e-01, -5.44834790e-01,  2.85788188e+00,
         -2.19098508e+00,  4.60275312e+00, -8.65932060e-01,
          0.00000000e+00,  4.60863886e-01,  3.86827298e+00,
          3.06025874e+00,  2.76800587e+00,  5.07200064e+00,
         -1.78508568e+00,  3.47097688e+00,  1.35991717e+00,
         -1.01210059e+00,  5.30044932e+00,  8.40281322e-01,
         -6.63779580e-01,  1.89452782e+00]]

testC :: LA.Vector R
testC = LA.fromList [ -83.9752375 , -118.424966  ,  -46.176225  ,  -44.1028972 ,
         -43.8399505 , -244.68241   ,  -67.069879  , -121.004792  ,
        -107.077591  ,  -20.0569717 ,  -82.6716474 ,  -73.6981641 ,
        -119.077927  ,   -9.97899905,  -41.7855714 ,  -56.4305485 ,
         -97.9298325 ,  -73.3123735 ,  -51.2169115 , -104.712497  ,
         -75.722393  ,  -83.2776042 ,  -34.5447282 ,   -6.88224218,
         -52.5345935 , -137.090858  ,  -35.8500455 , -292.975139  ,
         -73.2142654 , -238.516364  , -140.177043  ,  -31.0575237 ,
        -103.671827  ,  -27.3962    ,  -46.8354687 ,  -75.2934877 ,
        -102.353034  ,  -29.844633  ,  -57.4895948 ,  -54.1848221 ,
         -52.8721653 ]

testB :: LA.Vector R
testB = LA.fromList [221.05926829, 169.5812565 ,  94.15852609, 163.24788766,
          80.63507598, 216.09500765, 138.62729537,  68.67400288,
         119.05053092, 143.89199893, 157.45733783,  66.5613261 ,
         101.26280664, 118.68390784, 167.11478687, 258.66914518,
         235.47902565, 149.35884445, 159.98860853, 153.67098065,
         156.57284101, 163.14217634, 165.19741617,  97.908423  ]

testTab :: LA.Matrix R
testTab = toTableau (-testC) testMat testB

testIntMask :: LA.Vector Bool
testIntMask = LA.fromList $ replicate 41 True

xTree :: BranchProblem
xTree = searchBBTreeMax $ constructParBranchAndCut obj testTab testIntMask (-testC)
main :: IO()
main = print xTree