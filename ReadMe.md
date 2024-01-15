Here, we provide some supplementary files for the paper[1]. 

- The 'Movie 1' shows the simulated deformation process of different samples. To download the video, one can click the item, then click "View raw". 
- Each 'Input files and UMAT' folder contains ABAQUS input files and the corresponding FORTRAN files (UMAT subroutines). The results presented in the paper can be duplicated by these files. 

[1] Li, Z., Wang, J., Hossain, M., & Kadapa, C. (2023). A general theoretical scheme for shape-programming of incompressible hyperelastic shells through differential growth. *International Journal of Solids and Structures*, *265–266*, 112128. https://doi.org/10.1016/j.ijsolstr.2023.112128

# Simulation(needs expertise)

Each folder also contains necessary files to run the simulation using ABAQUS. To run the simulation, please make sure 

- Visual Studio and Parallel Studio are linked to ABAQUS, such that it can be used for subroutine development. [Here](https://www.researchgate.net/publication/349991987_Linking_ABAQUS_20192020_and_Intel_oneAPI_Base_Toolkit_FORTRAN_Compiler) is a linking guide from the internet. 

- Submit the job through ABAQUS COMMAND window, for example

  ```fortran
  C:
  cd C:\$YourPath$\Input files and UMAT\Example1
  
  abaqus job=SweetMelon user=SweetMelon.for 
  ```

  

# Simulated results

### Example 1: surfaces of revolution

![Example 1](https://github.com/Jeff97/electro-growth-UEL-subroutine/blob/main/Example1.jpg)

### Example 2: the helical surface (Cereus Forbesii Spiralis)

![Example 2](https://github.com/Jeff97/electro-growth-UEL-subroutine/blob/main/MappingHelical.jpg)

### Example 3: helical rod (tendril of pumpkin)

![Example 3](https://github.com/Jeff97/electro-growth-UEL-subroutine/blob/main/MappingRod.jpg)
