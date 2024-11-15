#include <bits/stdc++.h>
#include "mpi.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void save_image(const char *filename, int ***answer, int width, int height)
{
    unsigned char *image_data = (unsigned char *)malloc(width * height * 3);
    for (int x = 0; x < width; ++x)
    {
        for (int y = 0; y < height; ++y)
        {
            image_data[(x * height + y) * 3 + 0] = (unsigned char)(answer[x][y][0] & 0xFF); // Red
            image_data[(x * height + y) * 3 + 1] = (unsigned char)(answer[x][y][1] & 0xFF); // Green
            image_data[(x * height + y) * 3 + 2] = (unsigned char)(answer[x][y][2] & 0xFF); // Blue
        }
    }
    stbi_write_png(filename, width, height, 3, image_data, width * 3);
    free(image_data);
}

float find_opacity(float value, float *voxel_values, float *opacity, int n)
{
    if (value <= voxel_values[0])
    {
        return opacity[0];
    }
    else if (value >= voxel_values[n - 1])
    {
        return opacity[n - 1];
    }
    for (int i = 0; i < n - 1; i++)
    {
        if (value >= voxel_values[i] && value <= voxel_values[i + 1])
        {
            float x = (value - voxel_values[i]) / (voxel_values[i + 1] - voxel_values[i]);
            return (x * opacity[i + 1] + (1 - x) * opacity[i]);
        }
    }
    return 0.0;
}

int find_color(int n, float value, float *values, float *colors)
{
    if (value <= values[0])
    {
        return (int)(255 * colors[0]);
    }
    else if (value >= values[n - 1])
    {
        return (int)(255 * colors[n - 1]);
    }
    for (int i = 0; i < n - 1; i++)
    {
        if (value >= values[i] && value <= values[i + 1])
        {
            float x = (value - values[i]) / (values[i + 1] - values[i]);
            return (int)(255 * (x * colors[i + 1] + (1 - x) * colors[i]));
        }
    }
    return 0;
}

int minimum(int a, int b)
{
    return (a < b) ? a : b;
}

double maximum(double a, double b)
{
    return (a > b) ? a : b;
}

int main(int argc, char *argv[])
{
    int numtasks, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double start_time = MPI_Wtime();

    char *file_name = argv[1];
    int split_x = atoi(argv[2]), split_y = atoi(argv[3]), split_z = atoi(argv[4]);
    float stepsize = 0.5;
    char *opacity_file_name;
    char *color_file_name;

    if (argc == 8)
    {
        stepsize = atof(argv[5]);
        opacity_file_name = argv[6];
        color_file_name = argv[7];
    }
    else
    {
        opacity_file_name = argv[5];
        color_file_name = argv[6];
    }

    FILE *color_file = fopen(color_file_name, "r");
    float values[100], red[100], green[100], blue[100];
    char comma;
    int color_count = 0;
    while (fscanf(color_file, "%f%c", &values[color_count], &comma) == 2)
    {
        fscanf(color_file, "%f%c", &red[color_count], &comma);
        fscanf(color_file, "%f%c", &green[color_count], &comma);
        fscanf(color_file, "%f%c", &blue[color_count], &comma);
        color_count++;
    }
    fclose(color_file);

    FILE *opacity_file = fopen(opacity_file_name, "r");
    float voxel_values[10], opacity_values[10];
    int op_count = 0;
    while (fscanf(opacity_file, "%f%c", &voxel_values[op_count], &comma) == 2)
    {
        fscanf(opacity_file, "%f%c", &opacity_values[op_count], &comma);
        op_count++;
    }
    fclose(opacity_file);

    int x_dim, y_dim, z_dim;
    sscanf(file_name, "Isabel_%dx%dx%d_float32.raw", &x_dim, &y_dim, &z_dim);
    int x_min = 0, x_max = x_dim, y_min = 0, y_max = y_dim, z_min = 0, z_max = z_dim;

    int x_length = (x_max - x_min) / split_x;
    int y_length = (y_max - y_min) / split_y;
    int z_length = (z_max - z_min) / split_z;

    if (myrank == 0)
    {
        double data_read_start = MPI_Wtime();
        FILE *file = fopen(file_name, "rb");

        // Allocating space for data
        float ***data = (float ***)malloc(x_dim * sizeof(float **));
        for (int i = 0; i < x_dim; i++)
        {
            data[i] = (float **)malloc(y_dim * sizeof(float *));
            for (int j = 0; j < y_dim; j++)
            {
                data[i][j] = (float *)malloc(z_dim * sizeof(float));
            }
        }

        // Reading data from file
        for (int k = 0; k < z_dim; k++)
        {
            for (int j = 0; j < y_dim; j++)
            {
                for (int i = 0; i < x_dim; i++)
                {
                    fread(&data[i][j][k], sizeof(float), 1, file);
                }
            }
        }
        double data_read_end = MPI_Wtime();
        fclose(file);

        // Distributing the domain
        double domain_decomposition_start = MPI_Wtime();
        for (int i = 1; i < numtasks; i++)
        {
            int start_x = x_min + (i % split_x) * x_length;
            int start_y = y_min + ((i / split_x) % split_y) * y_length;
            int start_z = z_min + (i / (split_x * split_y)) * z_length;
            int end_x = (start_x + x_length > x_max) ? x_max : start_x + x_length - 1;
            int end_y = (start_y + y_length > y_max) ? y_max : start_y + y_length - 1;
            int end_z = (start_z + z_length > z_max) ? z_max : start_z + z_length - 1;

            int data_size = (end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1);
            float *data_to_send = (float *)malloc(data_size * sizeof(float));

            int index = 0;
            for (int x = start_x; x <= end_x; x++)
            {
                for (int y = start_y; y <= end_y; y++)
                {
                    for (int z = start_z; z <= end_z; z++)
                    {
                        data_to_send[index++] = data[x][y][z];
                    }
                }
            }
            MPI_Send(data_to_send, data_size, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            free(data_to_send);
        }
        double domain_decomposition_end = MPI_Wtime();
        double domain_decomposition_time = domain_decomposition_end - domain_decomposition_start;

        int ***answer = (int ***)calloc((x_max - x_min), sizeof(int **));
        for (int i = 0; i < (x_max - x_min); i++)
        {
            answer[i] = (int **)calloc((y_max - y_min), sizeof(int *));
            for (int j = 0; j < (y_max - y_min); j++)
            {
                answer[i][j] = (int *)calloc(3, sizeof(int));
            }
        }

        float **opacity = (float **)calloc((x_max - x_min), sizeof(float *));
        for (int i = 0; i < (x_max - x_min); i++)
        {
            opacity[i] = (float *)calloc((y_max - y_min), sizeof(float));
        }

        double ray_casting_start = MPI_Wtime();
        // Processing in process 0
        int start_x = x_min, start_y = y_min, start_z = z_min;
        int end_x = (start_x + x_length > x_max) ? x_max : start_x + x_length - 1;
        int end_y = (start_y + y_length > y_max) ? y_max : start_y + y_length - 1;
        int end_z = (start_z + z_length > z_max) ? z_max : start_z + z_length - 1;

        for (int i = start_x; i <= end_x; i++)
        {
            for (int j = start_y; j <= end_y; j++)
            {
                float face_value = data[i][j][0];
                opacity[i][j] = 0;
                answer[i][j][0] = 0;
                answer[i][j][1] = 0;
                answer[i][j][2] = 0;

                for (float k = z_length - 1; k >= 0; k -= stepsize)
                {
                    int k_floor = (int)k;
                    int k_ceil = k_floor + 1;
                    float t = k - k_floor;

                    float interpolated_value = (1 - t) * data[i][j][k_floor] + t * data[i][j][k_ceil];
                    float curr_opacity = find_opacity(interpolated_value, voxel_values, opacity_values, op_count);
                    int curr_color = find_color(color_count, interpolated_value, values, red);
                    answer[i][j][0] = answer[i][j][0] * (1 - opacity[i][j]) + curr_color * curr_opacity;
                    curr_color = find_color(color_count, interpolated_value, values, green);
                    answer[i][j][1] = answer[i][j][1] * (1 - opacity[i][j]) + curr_color * curr_opacity;
                    curr_color = find_color(color_count, interpolated_value, values, blue);
                    answer[i][j][2] = answer[i][j][2] * (1 - opacity[i][j]) + curr_color * curr_opacity;
                    opacity[i][j] = curr_opacity;
                }
                answer[i][j][0] = (int)answer[i][j][0];
                answer[i][j][1] = (int)answer[i][j][1];
                answer[i][j][2] = (int)answer[i][j][2];
            }
        }
        double ray_casting_end = MPI_Wtime();
        double computation_time = ray_casting_end - ray_casting_start;

        for (int i = 0; i < x_dim; i++)
        {
            free(data[i]);
        }
        free(data);

        double total_composition_time = 0.0;
        // Collecting from other processes
        for (int i = 1; i < numtasks; i++)
        {
            double composition_start = MPI_Wtime();
            int start_x = x_min + (i % split_x) * x_length;
            int start_y = y_min + ((i / split_x) % split_y) * y_length;
            int start_z = z_min + (i / (split_x * split_y)) * z_length;
            int end_x = (start_x + x_length > x_max) ? x_max : start_x + x_length - 1;
            int end_y = (start_y + y_length > y_max) ? y_max : start_y + y_length - 1;
            int end_z = (start_z + z_length > z_max) ? z_max : start_z + z_length - 1;

            int *color_receieving = (int *)malloc(3 * (end_x - start_x + 1) * (end_y - start_y + 1) * sizeof(int));
            float *opacity_receieving = (float *)malloc((end_x - start_x + 1) * (end_y - start_y + 1) * sizeof(float));
            MPI_Recv(color_receieving, 3 * (end_x - start_x + 1) * (end_y - start_y + 1), MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(opacity_receieving, (end_x - start_x + 1) * (end_y - start_y + 1), MPI_FLOAT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int ptr_i = start_x, ptr_j = start_y;
            for (int d = 0; d < (end_x - start_x + 1) * (end_y - start_y + 1); d++)
            {
                if (start_z == 0)
                {
                    answer[ptr_i][ptr_j][0] = color_receieving[3 * d];
                    answer[ptr_i][ptr_j][1] = color_receieving[3 * d + 1];
                    answer[ptr_i][ptr_j][2] = color_receieving[3 * d + 2];
                    opacity[ptr_i][ptr_j] = opacity_receieving[d];
                    ptr_j = ptr_j + 1;
                    if (ptr_j == end_y + 1)
                    {
                        ptr_j = start_y;
                        ptr_i += 1;
                    }
                }
                else
                {
                    answer[ptr_i][ptr_j][0] = (int)(answer[ptr_i][ptr_j][0] + (1 - opacity[ptr_i][ptr_j]) * opacity_receieving[d] * color_receieving[3 * d]);
                    answer[ptr_i][ptr_j][1] = (int)(answer[ptr_i][ptr_j][1] + (1 - opacity[ptr_i][ptr_j]) * opacity_receieving[d] * color_receieving[3 * d + 1]);
                    answer[ptr_i][ptr_j][2] = (int)(answer[ptr_i][ptr_j][2] + (1 - opacity[ptr_i][ptr_j]) * opacity_receieving[d] * color_receieving[3 * d + 2]);
                    opacity[ptr_i][ptr_j] = opacity[ptr_i][ptr_j] + (1 - opacity[ptr_i][ptr_j]) * opacity_receieving[d];
                    answer[ptr_i][ptr_j][0] = minimum(answer[ptr_i][ptr_j][0], 255);
                    answer[ptr_i][ptr_j][1] = minimum(answer[ptr_i][ptr_j][1], 255);
                    answer[ptr_i][ptr_j][2] = minimum(answer[ptr_i][ptr_j][2], 255);
                    ptr_j = ptr_j + 1;
                    if (ptr_j == end_y + 1)
                    {
                        ptr_j = start_y;
                        ptr_i += 1;
                    }
                }
            }
            double composition_end = MPI_Wtime();
            total_composition_time += composition_end - composition_start;
            free(color_receieving);
            free(opacity_receieving);
        }

        for (int i = 0; i < (end_x - start_x); i++)
        {
            free(opacity[i]);
        }
        free(opacity);

        char final_image_name[] = "%d_%d_%d.png";
        sprintf(final_image_name, final_image_name, split_x, split_y, split_z);
        save_image(final_image_name, answer, x_max - x_min, y_max - y_min);

        for (int i = 0; i < y_max - y_min; i++)
        {
            free(answer[i]);
        }
        free(answer);

        double end_time = MPI_Wtime();

        double comp_time;
        for (int i = 1; i < numtasks; i++)
        {
            MPI_Recv(&comp_time, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            computation_time = maximum(computation_time, comp_time);
        }
        printf("Communication time: %f\n", domain_decomposition_time + total_composition_time);
        printf("Computation time: %f\n", computation_time);
        printf("Total time: %f\n", end_time - start_time);
    }

    else
    {
        int start_x = x_min + (myrank % split_x) * x_length;
        int start_y = y_min + ((myrank / split_x) % split_y) * y_length;
        int start_z = z_min + (myrank / (split_x * split_y)) * z_length;
        int end_x = (start_x + x_length > x_max) ? x_max : start_x + x_length - 1;
        int end_y = (start_y + y_length > y_max) ? y_max : start_y + y_length - 1;
        int end_z = (start_z + z_length > z_max) ? z_max : start_z + z_length - 1;

        int data_size = (end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1);
        int *color_vector = (int *)malloc(3 * (end_x - start_x + 1) * (end_y - start_y + 1) * sizeof(int));
        float *opacity_vector = (float *)malloc((end_x - start_x + 1) * (end_y - start_y + 1) * sizeof(float));

        float *data_to_receive = (float *)malloc(data_size * sizeof(float));
        MPI_Recv(data_to_receive, data_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double ray_casting_start = MPI_Wtime();
        for (int d = 0; d < (end_x - start_x + 1) * (end_y - start_y + 1); d++)
        {
            float face_value = data_to_receive[d * (end_z - start_z + 1)];
            color_vector[3 * d] = 0;
            color_vector[3 * d + 1] = 0;
            color_vector[3 * d + 2] = 0;
            opacity_vector[d] = 0;

            for (float k = (end_z - start_z); k >= 0; k -= stepsize)
            {
                int k_floor = (int)k;
                int k_ceil = k_floor + 1;
                float t = k - k_floor;

                float interpolated_value = (1 - t) * data_to_receive[d * (end_z - start_z + 1) + k_floor] + t * data_to_receive[d * (end_z - start_z + 1) + k_ceil];
                int curr_color = find_color(color_count, interpolated_value, values, red);
                float curr_opacity = find_opacity(interpolated_value, voxel_values, opacity_values, op_count);

                color_vector[3 * d] = color_vector[3 * d] * (1 - opacity_vector[d]) + curr_color * curr_opacity;
                curr_color = find_color(color_count, interpolated_value, values, green);
                color_vector[3 * d + 1] = color_vector[3 * d + 1] * (1 - opacity_vector[d]) + curr_color * curr_opacity;
                curr_color = find_color(color_count, interpolated_value, values, blue);
                color_vector[3 * d + 2] = color_vector[3 * d + 2] * (1 - opacity_vector[d]) + curr_color * curr_opacity;
                opacity_vector[d] = curr_opacity;
            }
            color_vector[3 * d] = (int)color_vector[3 * d];
            color_vector[3 * d + 1] = (int)color_vector[3 * d + 1];
            color_vector[3 * d + 2] = (int)color_vector[3 * d + 2];
        }

        MPI_Send(color_vector, 3 * (end_x - start_x + 1) * (end_y - start_y + 1), MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(opacity_vector, (end_x - start_x + 1) * (end_y - start_y + 1), MPI_FLOAT, 0, 2, MPI_COMM_WORLD);

        double ray_casting_end = MPI_Wtime();
        double comp_time = ray_casting_end - ray_casting_start;
        MPI_Send(&comp_time, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

        free(opacity_vector);
        free(data_to_receive);
    }

    MPI_Finalize();
    return 0;
}
