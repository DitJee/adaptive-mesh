# adaptive mesh (staging)
The branch is called `staging` which will be used as a 'sandbox' for the project. This is where all changes from other branches except `main` will merge to.

# How to contribute to staging
To commit changes to `staging`, you must create a new branch based on `staging` where the changes will be made. After the changes are finished, you have to create a pull request to from `base`: `[your new branch]` to `master`: `staging` and wait for the code to be reviewed and approved.

# adaptive mesh
 This respository is about converting adaptive mesh code written in Fortran into Matlab. The theoretical of this project is based on Finite volume method in chapter 10 of Compressible Flow written by Prof. Pramote Dechaumphai, Ph.D.

 - Finite volume code, `hiflow`, in chapter 10 was initially written in Fortran. Then, the professor converted it into Matlab code since it is easier to develop and debug. Moreover, by developing in Matlab, the mesh can be plotted to verify the model and instantly display the results.
 - The code is now used to teach 3rd year students, and the data files which are mainly meshes have to be prepared before hand which is very tedious task.
 - To generate these mesh files, the code has to be developed for every single mesh files to many generate small cells.
 - The adaptive mesh is to way to automatically generate smaller cells where the shock waves occured to capture the shock wave and bigger cells in other areas to reduce the computation time.
 - This remeshing code is called `remesh` which is developed in Fortran.
 - The professor have grad student run the program using Fortran compiler. The result is documented in `"\adaptive-mesh\Info files\Comments on Remesh & Hiflow Codes.docx"`. It was stated that the problem with Fortran is that it is not suitable for displaying results.
 - So, the purpose of this project is to make future students be able to conveniently generate adaptive mesh and plot the result and solutions.

MIND IS HERE