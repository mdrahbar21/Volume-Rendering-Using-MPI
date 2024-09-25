#include <bits/stdc++.h>
#include "mpi.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void save_image(const char *filename, int ***answer, int width, int height)
{
    unsigned char *image_data = (unsigned char *)malloc(width * height * 3);
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            image_data[(y * width + x) * 3 + 0] = (unsigned char)(answer[y][x][0] & 0xFF); // Red
            image_data[(y * width + x) * 3 + 1] = (unsigned char)(answer[y][x][1] & 0xFF); // Green
            image_data[(y * width + x) * 3 + 2] = (unsigned char)(answer[y][x][2] & 0xFF); // Blue
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

int main(int argc, char *argv[])
{
    int numtasks, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    char *file_name = argv[1];
    int partition_dimension = atoi(argv[2]);
    float stepsize = 0.5;
    int x_min, x_max, y_min, y_max, x_dim, y_dim, z_dim;

    if (argc == 8)
    {
        stepsize = atof(argv[3]);
        x_min = atoi(argv[4]);
        x_max = atoi(argv[5]);
        y_min = atoi(argv[6]);
        y_max = atoi(argv[7]);
    }
    else
    {
        x_min = atoi(argv[3]);
        x_max = atoi(argv[4]);
        y_min = atoi(argv[5]);
        y_max = atoi(argv[6]);
    }
    sscanf(file_name, "Isabel_%dx%dx%d_float32.raw", &x_dim, &y_dim, &z_dim);

    float opacity_threshold = 0.97;

    int split_x = numtasks - 1, split_y = 1;
    if (partition_dimension == 2)
    {
        for (int i = 2; i < numtasks - 1; i++)
        {
            if ((numtasks - 1) % i == 0)
            {
                split_y = i;
                split_x = (numtasks - 1) / i;
                break;
            }
        }
    }

    int x_length = (x_max - x_min) / split_x;
    int y_length = (y_max - y_min) / split_y;

    if (myrank == 0)
    {
        FILE *file = fopen(file_name, "rb");

        float ***data = (float ***)malloc(x_dim * sizeof(float **));
        for (int i = 0; i < x_dim; i++)
        {
            data[i] = (float **)malloc(y_dim * sizeof(float *));
            for (int j = 0; j < y_dim; j++)
            {
                data[i][j] = (float *)malloc(z_dim * sizeof(float));
            }
        }

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

        for (int i = 0; i < numtasks - 1; i++)
        {
            int start_x = x_min + (i % split_x) * x_length;
            int start_y = y_min + (i / split_x) * y_length;
            int end_x = (i % split_x == split_x - 1) ? x_max : start_x + x_length;
            int end_y = (i / split_x == split_y - 1) ? y_max : start_y + y_length;

            float *data_to_send = (float *)malloc(z_dim * (end_x - start_x) * (end_y - start_y) * sizeof(float));
            int ptr_i = start_x, ptr_j = start_y, ptr_z = 0;
            for (int d = 0; d < z_dim * (end_x - start_x) * (end_y - start_y); d++)
            {
                data_to_send[d] = data[ptr_i][ptr_j][ptr_z];
                ptr_z += 1;

                if (ptr_z == z_dim)
                {
                    ptr_i = ptr_i + 1;
                    ptr_z = 0;
                    if (ptr_i == end_x)
                    {
                        ptr_i = start_x;
                        ptr_j += 1;
                    }
                }
            }
            MPI_Send(data_to_send, z_dim * (end_x - start_x) * (end_y - start_y), MPI_FLOAT, i + 1, 0, MPI_COMM_WORLD);
            free(data_to_send);
        }
        for (int i = 0; i < x_dim; i++)
        {
            free(data[i]);
        }
        free(data);
        fclose(file);

        int ***answer = (int ***)malloc((y_max - y_min) * sizeof(int **));
        for (int i = 0; i < (y_max - y_min); i++)
        {
            answer[i] = (int **)malloc((x_max - x_min) * sizeof(int *));
            for (int j = 0; j < (x_max - x_min); j++)
            {
                answer[i][j] = (int *)malloc(3 * sizeof(int));
            }
        }

        int total_opaque = 0;
        for (int i = 0; i < numtasks - 1; i++)
        {
            int start_x = (i % split_x) * x_length;
            int start_y = (i / split_x) * y_length;
            int end_x = (i % split_x == split_x - 1) ? (x_max - x_min) : start_x + x_length;
            int end_y = (i / split_x == split_y - 1) ? (y_max - y_min) : start_y + y_length;

            int *data_receieving = (int *)malloc(3 * (end_x - start_x) * (end_y - start_y) * sizeof(int));
            MPI_Recv(data_receieving, 3 * (end_x - start_x) * (end_y - start_y), MPI_INT, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int opaque_here;
            MPI_Recv(&opaque_here, 1, MPI_INT, i + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_opaque += opaque_here;
            int ptr_i = start_x, ptr_j = start_y;
            for (int d = 0; d < (end_x - start_x) * (end_y - start_y); d++)
            {
                answer[ptr_j][ptr_i][0] = data_receieving[3 * d];
                answer[ptr_j][ptr_i][1] = data_receieving[3 * d + 1];
                answer[ptr_j][ptr_i][2] = data_receieving[3 * d + 2];
                ptr_i = ptr_i + 1;
                if (ptr_i == end_x)
                {
                    ptr_i = start_x;
                    ptr_j += 1;
                }
            }
            free(data_receieving);
        }

        save_image("output.png", answer, x_max - x_min, y_max - y_min);

        for (int i = 0; i < y_max - y_min; i++)
        {
            free(answer[i]);
        }
        free(answer);
        float percent_opaque = total_opaque;
        percent_opaque /= ((x_max - x_min) * (y_max - y_min));
        printf("Total percentage of terminated rays is %.2f", 100 * percent_opaque);
    }

    else
    {
        // printf("Hello from task number %d\n", myrank);
        FILE *color_file = fopen("color_TF.txt", "r");
        float values[100], red[100], green[100], blue[100];
        char comma;
        int count = 0;
        while (fscanf(color_file, "%f%c", &values[count], &comma) == 2)
        {
            fscanf(color_file, "%f%c", &red[count], &comma);
            fscanf(color_file, "%f%c", &green[count], &comma);
            fscanf(color_file, "%f%c", &blue[count], &comma);
            count++;
        }
        // printf("count is %d\n", count);
        fclose(color_file);

        FILE *opacity_file = fopen("opacity_TF.txt", "r");
        float voxel_values[10], opacity_values[10];
        int op_count = 0;
        while (fscanf(opacity_file, "%f%c", &voxel_values[op_count], &comma) == 2)
        {
            fscanf(opacity_file, "%f%c", &opacity_values[op_count], &comma);
            op_count++;
        }
        fclose(opacity_file);

        int start_x = x_min + ((myrank - 1) % split_x) * x_length;
        int start_y = y_min + ((myrank - 1) / split_x) * y_length;
        int end_x = ((myrank - 1) % split_x == split_x - 1) ? x_max : start_x + x_length;
        int end_y = ((myrank - 1) / split_x == split_y - 1) ? y_max : start_y + y_length;

        int *color = (int *)malloc(3 * (end_x - start_x) * (end_y - start_y) * sizeof(int));
        float *opacity = (float *)malloc((end_x - start_x) * (end_y - start_y) * sizeof(float));

        for (int d = 0; d < (end_x - start_x) * (end_y - start_y); d++)
        {
            color[3 * d] = 0;
            color[3 * d + 1] = 0;
            color[3 * d + 2] = 0;
            opacity[d] = 0.0;
        }
        int num_opaque = 0;

        float *data_to_receive = (float *)malloc(z_dim * (end_x - start_x) * (end_y - start_y) * sizeof(float));
        MPI_Recv(data_to_receive, z_dim * (end_x - start_x) * (end_y - start_y), MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int d = 0; d < (end_x - start_x) * (end_y - start_y); d++)
        {
            for (float k = 0.0; k < z_dim - 1; k += stepsize)
            {
                int k_floor = (int)k;
                int k_ceil = k_floor + 1;
                float t = k - k_floor;

                float interpolated_value = (1 - t) * data_to_receive[d * z_dim + k_floor] + t * data_to_receive[d * z_dim + k_ceil];
                int curr_color = find_color(count, interpolated_value, values, red);
                float curr_opacity = find_opacity(interpolated_value, voxel_values, opacity_values, op_count) * stepsize;
                color[3 * d] = (int)(color[3 * d] + (1 - opacity[d]) * curr_color * curr_opacity);
                curr_color = find_color(count, interpolated_value, values, green);
                color[3 * d + 1] = (int)(color[3 * d + 1] + (1 - opacity[d]) * curr_color * curr_opacity);
                curr_color = find_color(count, interpolated_value, values, blue);
                color[3 * d + 2] = (int)(color[3 * d + 2] + (1 - opacity[d]) * curr_color * curr_opacity);
                opacity[d] = opacity[d] + (1 - opacity[d]) * curr_opacity;

                if (opacity[d] >= opacity_threshold)
                {
                    num_opaque += 1;
                    break;
                }
            }
        }
        MPI_Send(color, 3 * (end_x - start_x) * (end_y - start_y), MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&num_opaque, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        free(opacity);
        free(data_to_receive);
    }

    MPI_Finalize();
    return 0;
}
