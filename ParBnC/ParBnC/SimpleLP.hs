{-# LANGUAGE RecordWildCards #-}
import Numeric.LinearProgramming
    ( (#),
      simplex,
      Bound((:==:)),
      Constraints(Sparse, Dense),
      Optimization(Maximize, Minimize),
      Solution, Bounds, Bound( (:&:) ))

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

