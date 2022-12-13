{-# LANGUAGE RecordWildCards #-}
import Numeric.LinearProgramming
    ( (#),
      simplex,
      Bound((:==:)),
      Constraints(Sparse, Dense),
      Optimization(Maximize, Minimize),
      Solution, Bounds, Bound( (:&:) ))
import Numeric.IEEE ( IEEE(infinity) )
import Data.List
import Numeric.LinearAlgebra as LA

prob :: Optimization
prob = Maximize [4, -3, 2, 0, 0]

constr1 :: Constraints
constr1 = Sparse [ [2#1, 1#2, 1#4] :==: 10
                 , [1#2, 5#3, 1#5] :==: 20
                 ]

-- simplex prob constr1 []

data MatConstraints = MatVec [[Double]] [Double] deriving Show

data VariableType = INTEGER | CONTINUOUS deriving Show

data ProblemType = LP | MIP deriving Show

data Problem = Problem { 
    problemType :: ProblemType,
    variables :: [String],
    objective :: Optimization,
    constraints :: MatConstraints,
    lowerBounds :: [Double],
    upperBounds :: [Double],
    variableTypes :: [VariableType]
}

data ObjectiveType = Maximization | Minimization deriving Show

parseConstraints :: MatConstraints -> Constraints
parseConstraints (MatVec aMatrix bVector) = Dense $ zipWith (:==:) aMatrix bVector

parseBounds :: [Double] -> [Double] -> Bounds
parseBounds lowerBounds upperBounds = zipWith3 f index lowerBounds upperBounds where
    n = length lowerBounds
    index = [1 .. n]
    f x y z = x :&: (y, z)

solveProblem :: Problem -> Solution
solveProblem Problem {..} = simplex objective constr bnd where
    constr = parseConstraints constraints
    bnd = parseBounds lowerBounds upperBounds

branchAndBoundSolve :: Problem -> Solution
branchAndBoundSolve x@Problem {problemType = LP, ..} = solveProblem x

-- for simple LP test


toTableau :: LA.Vector R -> LA.Matrix R -> LA.Vector R -> LA.Matrix R
toTableau costC matA constB = tab where 
    xb = LA.fromColumns $ LA.toColumns matA ++ [constB]
    z = mappend costC $ vector [0]
    tab = LA.fromRows $ LA.toRows xb ++ [z]

costCheck :: ObjectiveType -> (R -> Bool)
costCheck Maximization = (> 0)
costCheck Minimization = (< 0)

isImprovable :: ObjectiveType -> LA.Matrix R -> Bool
isImprovable obj tab = any (costCheck obj) $ LA.toList cost where 
    cost = subVector 0 (cols tab - 1) $ tab ! (rows tab - 1)

getPivotPosition :: ObjectiveType -> LA.Matrix R -> (Int, Int)
getPivotPosition obj tab = (row, column) where 
    z = tab ! (rows tab - 1)
    cost = subVector 0 (cols tab - 1) z
    column = head $ LA.find (costCheck obj) cost
    getElem rowEq
        | elem <= 0 = infinity::R
        | otherwise = (rowEq ! (LA.size rowEq - 1)) / elem
        where 
            elem = rowEq ! column
    restrictions = map getElem $ init (LA.toRows tab)
    Just row =  elemIndex (minimum restrictions) restrictions

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
updateTab obj subTab
    | not $ isImprovable obj subTab = subTab
    | otherwise = updateTab obj newTab where 
        pivotPos = getPivotPosition obj subTab
        newTab = pivotStep subTab pivotPos

simplexWithTab :: ObjectiveType -> LA.Vector R -> LA.Matrix R -> LA.Vector R -> LA.Vector R
simplexWithTab obj costC matA constB = solution where 
    tab = toTableau costC matA constB       
    lastTab = updateTab obj tab
    solution = getSolution lastTab

b :: Vector R
b = LA.fromList [2,4,4::R]
matA :: Matrix R
matA = LA.fromLists [[-1,1,1,0,0],[1,0,0,1,0],[0,1,0,0,1::R]]
c :: Vector R
c = fromList [1,1,0,0,0::R]


problemTest :: Problem
problemTest = Problem { 
    problemType = LP,
    variables = [" "],
    objective = Maximize [4, -3, 2, 0, 0],
    constraints = MatVec [[2,1,0,1,0],[0,1,5,0,1]] [10,20],
    lowerBounds = [0,0,0,0,0],
    upperBounds = [100,100,100,100,100],
    variableTypes = [CONTINUOUS| _ <- [0..5]]
}

