# CS677 Topics in Large Data Analysis and Visualization

# Assignment 1

### Group Members (Group 7)

1. Arindom Bora (210183)
2. Khush Khandelwal (210511)
3. Md Rahbar (210601)

### Instructions to run the code

1. Ensure both the <code>g7solution.cpp</code> and <code>stb_image_write.h</code> are in the same folder.
2. Compile the code using <code>mpicxx -o executable_file_name g7solution.cpp</code>
3. Execute the code giving appropriate arguments. For example,
   <code>mpirun -np 8 ./executable_file_name Isabel_1000x1000x200_float32.raw 2 0.25 0 500 0 999</code>

### Overview of the code

This program is designed to read volumetric image data, partition it across multiple tasks using MPI, process the image slices in parallel, and generate a final image output (PNG format). The program handles color mapping and opacity calculations on voxel data and saves the rendered image.

The code can be divided into two major components:

1. **Root Process (rank 0):** Responsible for reading the volumetric image data, distributing the workload to worker processes, collecting the processed results, and saving the final image.
2. **Worker Processes:** Each worker process computes its assigned section of the volumetric data using color and opacity mappings, then sends the processed image data back to the root.

#### Dependencies

1. **MPI:** Used for parallel processing
2. **stb_image_write:** A lightweight library to save the processed image in PNG format which can be downloaded from https://github.com/nothings/stb/blob/root/stb_image_write.h

#### Important functions

**1. Save 3D matrix of RGB values into PNG image**
A function **save_image** which takes inputs of the final image name, the color matrix, the width and height of the image size. It converts the integer RGB values into unsigned char and writes them to a PNG image using predefined functions from the header file **stb_image_write.h**
**2. Color and Opacity Mapping**
Given a float value of the current voxel, these functions run a linear search on the values present in the Transfer function and then perform linear interpolation to find the desried color/ opacity value.

#### Domain Decomposition

Based on the input for 1D or 2D partitioning, the code assign the number of splits to happen in each direction. If n is the total number of parallel processes, for 1D partitioning, number of splits in X-direction = n-1 (excluding root Process), while for 2D partitioning, the number of splits in the Y-direction is taken to be the least factor of n-1 greater than 1 (1 if such factor doesn't exist).

### Reading Data

The data is read using 3 for loops and data is filled in the X-direction first, followed by Y and then Z.

### Computation

**1. Send data over MPI**
In order to send data from the root process to the worker process, the data is first flattened to a 1D array where the values in Z-dimension are filled first, followed by Y and then X.
**2. Receive data in worker processes**
The flattened array is received in each of the worker processes where ray casting is done in the z-direction. For fractional step-size, the values are linearly interpolated and then color and opacity values are calculated. The opacity is multiplied with the step-size to prevent it from converging faster.
**3. Early Ray Termination**
A threshold of 0.97 is kept. If any ray crosses this opacity value, the ray is terminated and such early terminated rays are counted.
**4. Send Color and ray termination data from worker to root**
A flattend 2D array of RGB values are sent to the root process from each worker process where they are accumulated over a 3D array.
Also, the total percentage of terminated rays is calculated in the root process by summing over the number of terminated rays in the worker processes.
