# CS 677 Assignment 2: Parallel 3D Volume Rendering

### Group Members (Group 7)

1. Arindom Bora (210183)
2. Khush Khandelwal (210511)
3. Md Rahbar (210601)

---

### Instructions to run the code

1. Ensure the <code>solution.cpp</code>, <code>stb_image_write.h</code>, <code>color_TF.txt</code> and <code>opacity_TF.txt</code> files are in the same folder.
2. Compile the code using <code>mpicxx -o executable_file_name solution.cpp</code>
3. Execute the code giving appropriate arguments. For example,<br/>
   <code>mpirun -np 8 ./executable Isabel_1000x1000x200_float32.raw 2 2 2 0.5 opacity_TF.txt color_TF.txt</code><br/>or<br/>
   <code>mpirun-np 16 --host csews16:8,csews19:8 --mca btl_tcp_if_include eno1 --mca btl tcp,self ./executable Isabel_1000x1000x200_float32.raw 2 2 2 0.5 opacity_TF.txt color_TF.txt</code>

   Use --mca to run code on common network interface, if there are multiple network interfaces among the servers. In the above command, I am running on "eno1" interface.

---

#### Dependencies

1. **MPI:** Used for parallel processing
2. **stb_image_write:** A lightweight library to save the processed image in PNG format which can be downloaded from https://github.com/nothings/stb/blob/root/stb_image_write.h

---
