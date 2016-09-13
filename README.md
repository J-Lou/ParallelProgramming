# ParallelProgramming

Assignment overview
-----------------------------

The main goal of this assignment is to gain the experience of developing a parallel program in pthreads. Additionally, you will study how the parallel performance is affected by your parallelization strategy (including task decomposition, assignment, and synchronization management) and the multiprocessor platform.

Requirements for the main assignment
-----------

You are asked to develop a parallel Gaussian Elimination (with partial pivoting) program using pthreads. In the context of solving a system of linear equations, Gaussian Elimination is a classical method for reducing an equation matrix into an equivalent upper-diagonal matrix. The full solution requires an additional back-solver, which has a lower computation complexity and thus typically requires much less time than the Gaussian Elimination step. In this assignment, you only need to parallelize the Gaussian Elimination step.
A sequential version of Gaussian Elimination is provided to you. It is available at Gauss/seq . The program takes one parameter as the input matrix. Currently six matrices (five real-world sparse matrices and one artificially generated dense matrix) are available for your testing. Please read Gauss/README for details.

After developing and testing your program, you should measure the performance/speedup of your program with the input matrices on at least two different multiprocessor machines (for up to at least four processors). You should analyze the performance results. In addition, I expect that you have tried multiple ways of realizing parallel Gaussian Elimination (in terms of task decomposition, assignment, or synchronization management). You need to provide a comparison on at least two alternative approaches (pick a comparison that you have learned the most from) and analyze their performance.

We only care about the timing of the Gaussian Elimination step of your program. In your parallel program, please use barrier (see the sor code) to properly synchronize all processors at the beginning and end of Gaussian Elimination for timing. Code for using the high resolution timer on Intel processors is available at . To acquire accurate timing, you'll need to make sure that the machine has enough number of processors for your experiments.