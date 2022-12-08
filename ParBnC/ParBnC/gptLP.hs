import Data.List

-- function to find the pivot element in the given row and column
findPivot :: [[Int]] -> Int -> Int -> (Int, Int)
findPivot matrix row col = let
-- get the column from the given matrix
    column = map (!! col) matrix
-- find the minimum non-negative element in the column
    pivot = minimum $ filter (>=0) column
-- get the index of the pivot element in the column
    Just pivotIndex = elemIndex pivot column
    in
-- return the row and column indices of the pivot element
    (pivotIndex, col)

-- function to perform a row operation on the given matrix
performRowOperation :: [[Int]] -> (Int, Int) -> [[Int]]
performRowOperation matrix (row, col) = let
-- get the pivot element
    pivot = matrix !! row !! col
-- divide each element in the pivot row by the pivot element
    pivotRow = map (`div` pivot) $ matrix !! row
-- get the rest of the matrix
    rest = take row matrix ++ drop (row + 1) matrix
-- update the matrix by replacing the pivot row with the divided row
    updatedMatrix = take row matrix ++ [pivotRow] ++ drop row matrix
    in
    updatedMatrix

-- function to perform a column operation on the given matrix
performColumnOperation :: [[Int]] -> (Int, Int) -> [[Int]]
performColumnOperation matrix (row, col) = let
-- get the pivot element
    pivot = matrix !! row !! col
    -- update each row in the matrix by subtracting the pivot row multiplied by the corresponding element in the pivot column
    updatedMatrix = map (\x -> zipWith (-) x $ map (* pivot) (matrix !! row)) matrix
    in
    updatedMatrix

-- function to solve the given linear program
solveLinearProgram :: [[Int]] -> [Int] -> [Int] -> (Int, [Int])
solveLinearProgram matrix objectiveCoeffs constraints = let
-- get the number of variables in the linear program
    numVariables = length objectiveCoeffs
    -- get the number of constraints in the linear program
    numConstraints = length constraints
    -- create an initial tableau by appending the objective function coefficients and the constraints to the matrix
    tableau = matrix ++ [objectiveCoeffs] ++ map ((\ x -> x ++ [1]) . zipWith (*) constraints) matrix
    -- find the pivot element in the first row and last column
    pivot = findPivot tableau 0 (numVariables + numConstraints)
    -- perform row and column operations on the tableau until the optimal solution is reached
    optimalTableau = iterate (\m -> performColumnOperation (performRowOperation m pivot) pivot) tableau !! numVariables
    -- extract the optimal solution from the final tableau
    optimalSolution = map last $ take numConstraints optimalTableau
    in
    -- return the optimal value and the optimal solution
    (last $ last optimalTableau, optimalSolution)

mat :: [[Int]]
mat = [[1,0,0], [2,1,-1]]
obj :: [Int]
obj = [1,1,1]
constraint :: [Int]
constraint = [4,3]