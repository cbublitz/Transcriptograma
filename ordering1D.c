#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct graph_ptr{
   short cur;
   short ini;
} graph_ptr;

typedef struct args{
   float alpha;
   char *file_name;
   float cooling_factor;
   short cooling_interval;
   float percentual_energy;
   unsigned short mc_steps;
   double seed;
} args;

args args_parser(int argc, char *argv[]){
   args options;
   short x;

   options.alpha = 1.0;
   options.mc_steps = 100;
   options.file_name = NULL;
   options.seed = time(NULL);
   options.cooling_factor = 0.5;
   options.cooling_interval = 100;
   options.percentual_energy = 0.0001;

   /* input parameters */
   printf("Parameters list:\n");
   printf("\tf=Association file\n");
   printf("\ti=Number of isothermal steps\n");
   printf("\tm=Number of Monte Carlo steps\n");
   printf("\tc=Cooling factor\n");
   printf("\ta=Alpha value\n");
   printf("\tp=Percentual energy for initial temperature\n");
   printf("\ts=Random seed\n\n");

   for(x=1;x<argc;x++){
      switch(argv[x][0]){
         case 'f':
            options.file_name = &argv[x][2];
            break;
         case 'i':
            options.cooling_interval = atoi(&argv[x][2]);
            break;
         case 'm':
            options.mc_steps = atoi(&argv[x][2]);
            break;
         case 'c':
            options.cooling_factor = atof(&argv[x][2]);
            break;
         case 'a':
            options.alpha = atof(&argv[x][2]);
            break;
         case 'p':
            options.percentual_energy = atof(&argv[x][2]);
            break;
         case 's':
            options.seed = atoi(&argv[x][2]);
            break;
      }
   }

   return options;
}

void swap(graph_ptr *order, short a, short b){
   short aux_swap;

   aux_swap = order[order[a].cur].ini;
   order[order[a].cur].ini = order[order[b].cur].ini;
   order[order[b].cur].ini = aux_swap;

   aux_swap = order[a].cur;
   order[a].cur = order[b].cur;
   order[b].cur = aux_swap;
}

void get_size_column(char *matrix[], short *size_column, unsigned short n_nodes){
   short col, row;

   for(row=(n_nodes-1);row>=0;row--){
      size_column[row] = 0;
      for(col=0;col<n_nodes;col++){
         if(matrix[row][col] == 1){
            size_column[row]++;
         }
      }
   }
}

void compress_matrix(char *matrix[], short *ones_list[], unsigned short n_nodes){
   short col, row, aux_col;

   for(row=0;row<n_nodes;row++){
      aux_col = 0;
      for(col=0;col<n_nodes;col++){
         if(matrix[row][col] == 1){
            ones_list[row][aux_col] = col;
            aux_col++;
         }
      }
   }
}

double getPartEnergy(char *matrix[], graph_ptr *order, short a, short b, short *size_column, short *ones_list[], unsigned short n_nodes, float alpha){
   short index_a = a, index_b = b, init_a, init_b, row, col, index_cnt;
   double col_a_energy = 0, col_b_energy = 0, row_a_energy = 0, row_b_energy = 0, part_energy, abs;
   unsigned char index_a_cnt, index_b_cnt, neighbors;
   graph_ptr pos_i, pos_j;

   if(a > b){
      index_a = b;
      index_b = a;
   }

   init_a = index_a;
   init_b = index_b;

   if((index_a+1) == index_b){
      if(index_a == 0){
         index_a_cnt = 2;
      }else{
         index_a_cnt = 3;
         index_a = index_a-1;
      }

      if(index_b == (n_nodes-1)){
         index_b_cnt = 0;
      }else{
         index_b_cnt = 1;
         index_b = index_b+1;
      }

   }else if((index_a+1) == (index_b-1)){
      if(index_a == 0){
         index_a_cnt = 2;
      }else{
         index_a_cnt = 3;
         index_a = index_a-1;
      }

      if(index_b == (n_nodes-1))
         index_b_cnt = 1;
      else
         index_b_cnt = 2;

   }else{
      if(index_a == 0){
         index_a_cnt = 2;
      }else{
         index_a_cnt = 3;
         index_a = index_a-1;
      }

      if(index_b == (n_nodes-1))
         index_b_cnt = 2;
      else
         index_b_cnt = 3;

      index_b = index_b-1;
   }

   index_cnt = index_a + index_a_cnt;
   /* column a */
   for(col=index_a;col<index_cnt;col++){
      part_energy = 0;
      pos_i = order[col];
      for(row=0;row<size_column[pos_i.cur];row++){
         pos_j = order[ones_list[pos_i.cur][row]];
         if(col < pos_j.ini){
            neighbors = 0;
            if(pos_j.ini != 0)
               neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
            if(pos_j.ini != (n_nodes-1))
               neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
            if(col != 0)
               neighbors = neighbors + (1 - matrix[order[col-1].cur][order[pos_j.ini].cur]);
            if(col != (n_nodes-1))
               neighbors = neighbors + (1 - matrix[order[col+1].cur][order[pos_j.ini].cur]);

            abs = pow((pos_j.ini - col), alpha);
            part_energy = part_energy + (abs * neighbors);
         }
      }
      col_a_energy = col_a_energy + part_energy;
   }

   /* row a */
   for(row=(index_cnt-1);row>=index_a;row--){
      part_energy = 0;
      pos_i = order[row];
      for(col=0;col<size_column[pos_i.cur];col++){
         pos_j = order[ones_list[pos_i.cur][col]];
         if(row > pos_j.ini){
            if((pos_j.ini!= init_a-1) && (pos_j.ini != init_a) && (pos_j.ini != init_a+1) && (pos_j.ini != init_b-1) && (pos_j.ini != init_b) && (pos_j.ini != init_b+1)){
               neighbors = 0;
               if(pos_j.ini != 0)
                  neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
               if(pos_j.ini != (n_nodes-1))
                  neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
               if(row != 0)
                  neighbors = neighbors + (1 - matrix[order[row-1].cur][order[pos_j.ini].cur]);
               if(row != (n_nodes-1))
                  neighbors = neighbors + (1 - matrix[order[row+1].cur][order[pos_j.ini].cur]);

               abs = pow((row - pos_j.ini), alpha);
               part_energy = part_energy + (abs * neighbors);
            }
         }
      }
      row_a_energy = row_a_energy + part_energy;
   }

   if(index_b_cnt != 0){
      index_cnt = index_b + index_b_cnt;
      /* column b */
      for(col=index_b;col<index_cnt;col++){
         part_energy = 0;
         pos_i = order[col];
         for(row=0;row<size_column[pos_i.cur];row++){
            pos_j = order[ones_list[pos_i.cur][row]];
            if(col < pos_j.ini){
               neighbors = 0;
               if(pos_j.ini != 0)
                  neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
               if(pos_j.ini != (n_nodes-1))
                  neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
               if(col != 0)
                  neighbors = neighbors + (1 - matrix[order[col-1].cur][order[pos_j.ini].cur]);
               if(col != (n_nodes-1))
                  neighbors = neighbors + (1 - matrix[order[col+1].cur][order[pos_j.ini].cur]);

               abs = pow((pos_j.ini - col), alpha);
               part_energy = part_energy + (abs * neighbors);
            }
         }
         col_b_energy = col_b_energy + part_energy;
      }

      /* row b */
      for(row=(index_cnt-1);row>=index_b;row--){
         part_energy = 0;
         pos_i = order[row];
         for(col=0;col<size_column[pos_i.cur];col++){
            pos_j = order[ones_list[pos_i.cur][col]];
            if(row > pos_j.ini){
               if((pos_j.ini!= init_b-1) && (pos_j.ini != init_b) && (pos_j.ini != init_b+1) && (pos_j.ini != init_a-1) && (pos_j.ini != init_a) && (pos_j.ini != init_a+1)){
                  neighbors = 0;
                  if(pos_j.ini != 0)
                     neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
                  if(pos_j.ini != (n_nodes-1))
                     neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
                  if(row != 0)
                     neighbors = neighbors + (1 - matrix[order[row-1].cur][order[pos_j.ini].cur]);
                  if(row != (n_nodes-1))
                     neighbors = neighbors + (1 - matrix[order[row+1].cur][order[pos_j.ini].cur]);

                  abs = pow((row - pos_j.ini), alpha);
                  part_energy = part_energy + (abs * neighbors);
               }
            }
         }
         row_b_energy = row_b_energy + part_energy;
      }
   }
   return (col_a_energy + col_b_energy + row_a_energy + row_b_energy);
}

double getMatEnergy(char *matrix[], graph_ptr *order, short *size_column, short *ones_list[], unsigned short n_nodes, float alpha){
   double energy = 0, row_energy, abs;
   unsigned char neighbors;
   short row, col;
   graph_ptr pos_i, pos_j;

   for(row=(n_nodes-1);row>=0;row--){
      row_energy = 0;
      pos_i = order[row];
      for(col=0;col<size_column[pos_i.cur];col++){
         pos_j = order[ones_list[pos_i.cur][col]];
         if(row < pos_j.ini){
            neighbors = 0;
            if(pos_j.ini != 0)
               neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
            if(pos_j.ini != (n_nodes-1))
               neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
            if(row != 0)
               neighbors = neighbors + (1 - matrix[order[row-1].cur][order[pos_j.ini].cur]);
            if(row != (n_nodes-1))
               neighbors = neighbors + (1 - matrix[order[row+1].cur][order[pos_j.ini].cur]);

            abs = pow((pos_j.ini - row), alpha);
            row_energy = row_energy + (abs * neighbors);
         }
      }
      energy = energy + row_energy;
   }
   return energy;
}

void buildMat(FILE *inputfile, char ***matrix, char **nodes, unsigned short n_nodes){
   char line[80], *node_left, *node_right;
   /*short y, last_node = 0; unused*/
   short x, node_index_left, node_index_right;
   char **alloc_matrix = (char**) malloc(n_nodes * sizeof(char*));
   for(x=0;x<n_nodes;x++)
      alloc_matrix[x] = (char *) calloc(n_nodes,sizeof(char));

   while(fgets(line, 80, inputfile) != NULL){
      node_left = strtok(line," \t\n");
      node_right = strtok(NULL," \t\n");
      if((node_left != NULL) && (node_right != NULL)){
         node_index_left = -1;
         node_index_right = -1;
         for(x=0;x<n_nodes;x++){
            if(!strcmp(nodes[x],node_left))
               node_index_left = x;
            if(!strcmp(nodes[x],node_right))
               node_index_right = x;
            if((node_index_left>0) && (node_index_right>0))
               break;
         }
      }
      alloc_matrix[node_index_left][node_index_right] = 1;
      alloc_matrix[node_index_right][node_index_left] = 1;
   }
   *matrix = alloc_matrix;
}

unsigned short readProteinList(FILE *inputfile, char ***protein_list){
   char line[80], *node_left, *node_right;
   short x, last_node = 0, node_left_length, node_right_length;
   unsigned short left_found = 0, right_found = 0, n_nodes = 0;
   char **alloc_nodes = NULL;
   unsigned int n_edges = 0;

   while(fgets(line, 80, inputfile) != NULL){
      node_left = strtok(line," \t\n");
      node_left_length = strlen(node_left)+1;
      node_right = strtok(NULL," \t\n");
      node_right_length = strlen(node_right)+1;
      if((node_left != NULL) && (node_right != NULL)){
         if(n_nodes == 0){
            alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
            alloc_nodes[n_nodes] = (char *) calloc(node_left_length,sizeof(char));
            strcpy(alloc_nodes[n_nodes],node_left);
            left_found = 1;
            n_nodes++;
            alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
            alloc_nodes[n_nodes] = (char *) calloc(node_right_length,sizeof(char));
            strcpy(alloc_nodes[n_nodes],node_right);
            right_found = 1;
            n_nodes++;
         }else{
            if(!strcmp(alloc_nodes[last_node],node_left)){
               right_found = 0;
               for(x=0;x<n_nodes;x++){
                  if(!strcmp(alloc_nodes[x],node_right)){
                     right_found = 1;
                     break;
                  }else if(x == (n_nodes-1)){
                     alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
                     alloc_nodes[n_nodes] = (char *) calloc(node_right_length,sizeof(char));
                     strcpy(alloc_nodes[n_nodes],node_right);
                     right_found = 1;
                     n_nodes++;
                  }
               }
            }else{
               left_found = 0;
               right_found = 0;
               for(x=0;x<n_nodes;x++){
                  if(!strcmp(alloc_nodes[x],node_left) && !left_found)
                     left_found = 1;
                  if(!strcmp(alloc_nodes[x],node_right) && !right_found)
                     right_found = 1;
                  if(left_found && right_found)
                     break;
                  if(x == (n_nodes-1)){
                     if(!left_found){
                        alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
                        alloc_nodes[n_nodes] = (char *) calloc(node_left_length,sizeof(char));
                        strcpy(alloc_nodes[n_nodes],node_left);
                        last_node = n_nodes;
                        n_nodes++;
                     }
                     if(!right_found){
                        alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
                        alloc_nodes[n_nodes] = (char *) calloc(node_right_length,sizeof(char));
                        strcpy(alloc_nodes[n_nodes],node_right);
                        n_nodes++;
                     }
                  }
               }
            }
         }
         n_edges++;
      }
   }
   printf("\tProteins: %d\n",n_nodes);
   printf("\tInteractions: %d\n",n_edges/2);

   *protein_list = alloc_nodes;
   return n_nodes;
}

void rand_values(short *random_list, unsigned short n_nodes){
   short x, y, value;
   unsigned short repeat;

   for(x=0;x<n_nodes;x++){
      repeat = 1;
      random_list[x] = -1;
      do{
         value = (int)(rand()/(RAND_MAX + 1.0) * n_nodes);
         for(y=0;y<=x;y++){
            if((random_list[y]) == value){
               break;
            }else if(y == x){
               random_list[y] = value;
               repeat = 0;
            }
         }
      }while(repeat);
   }
}

void progressBar(unsigned short step, unsigned short monte_carlo_steps){
   float ratio;
   unsigned short incomplete, x;

   ratio = (float)step/monte_carlo_steps;
   incomplete = ratio * 100;

   printf("%3d%% [", (int)(ratio*100) );
   for(x=0;x<incomplete;x++)
      printf("=");
   for(x=incomplete;x<100;x++)
      printf(" ");
   printf("]\r");
}

/* ****************************************************************************************** */
/* -------------------------------------- MAIN PROGRAM -------------------------------------- */
/* ****************************************************************************************** */
int main (int argc, char *argv[]){
   args options;
   unsigned short step, n_nodes;
   graph_ptr *graph_order = NULL;
   short n, a, b, *random_list = NULL;
   double initial_energy, curr_energy, new_energy, delta_energy, curr_part_energy, new_part_energy, temperature, expo, swap_cnt = 0;
   long int path_chars;
   char **matrix = NULL, **protein_list = NULL, **nodes = NULL, *dir_name = NULL, *order_path = NULL, *dat_name = NULL;
   FILE *datfile, *order, *energy;

   options = args_parser(argc, argv);

   if(options.file_name != NULL){
      if(strrchr(options.file_name,'/') == NULL){
         path_chars = 0;
         dat_name = (char*) malloc((strlen(options.file_name) + 1)*sizeof(char));
         strncpy(dat_name,options.file_name,strlen(options.file_name));
         dat_name[strlen(options.file_name)] = '\0';
      }else{
         path_chars = strrchr(options.file_name,'/') - options.file_name;
         dat_name = (char*) malloc((strlen(strrchr(options.file_name,'/')))*sizeof(char));
         strncpy(dat_name,&options.file_name[path_chars + 1],strlen(strrchr(options.file_name,'/')) - 1);
         dat_name[strlen(strrchr(options.file_name,'/')) - 1] = '\0';
      }

      order_path = (char*) malloc((strlen(options.file_name) + 10)*sizeof(char));
      dir_name = (char*) malloc((path_chars + 1)*sizeof(char));
      strncpy(dir_name, options.file_name, path_chars);
      dir_name[path_chars] = '\0';

      datfile = fopen(options.file_name,"r");
      if(datfile == NULL){
         printf("ERROR! File %s does not exist.\n",options.file_name);
         exit(1);
      }else{
         printf("Reading file...\n");
         n_nodes = readProteinList(datfile, &protein_list);
         rewind(datfile);
      }
   }else{
      printf("ERROR! Invalid file name.\n");
      exit(1);
   }

   graph_order = (graph_ptr*) malloc(n_nodes*sizeof(graph_ptr));
   nodes = (char**) malloc(n_nodes*sizeof(char*));
   random_list = (short*) malloc(n_nodes*sizeof(short));

   srand(options.seed);
   printf("Applying random order...\n");
   rand_values(random_list, n_nodes);
   for(n=0;n<n_nodes;n++){
      nodes[random_list[n]] = protein_list[n];
      graph_order[n].cur = n;
      graph_order[n].ini = n;
   }
   free(random_list);
   free(protein_list);
   buildMat(datfile, &matrix, nodes, n_nodes);
   fclose(datfile);

   /* create a list with number of elements of each column */
   short size_column[n_nodes];
   get_size_column(matrix, size_column, n_nodes);

   /* create a list representation to compress the matrix */
   short **ones_list = (short**) malloc(n_nodes*sizeof(short*));
   for(n=0;n<n_nodes;n++){
      ones_list[n] = (short*) malloc(size_column[n]*sizeof(short));
   }
   compress_matrix(matrix, ones_list, n_nodes);
   initial_energy = getMatEnergy(matrix, graph_order, size_column, ones_list, n_nodes, options.alpha);
   printf("INITIAL ENERGY: %.f\n",initial_energy*2);
   curr_energy = initial_energy;
   temperature = curr_energy * options.percentual_energy;

   if(strrchr(dat_name,'.') == NULL)
      sprintf(order_path,"%s/energy_%s.dat",dir_name,dat_name);
   else
      sprintf(order_path,"%s/energy_%s",dir_name,dat_name);
   energy = fopen(order_path,"w+");
   if(energy == NULL){
      printf("ERROR! File %s does not exist.\n",order_path);
      exit(1);
   }
   fprintf(energy,"steps\tenergy\ttemperature\tswaps\n");
   fprintf(energy,"0\t%.f\t%.3lE\t%.lf\n",curr_energy*2, temperature*2, swap_cnt);

   printf("Ordering...\n");
   for(step=1;step<=options.mc_steps;step++){
      for(n=0;n<n_nodes;n++){
         /* get random nodes */
         a = (int)(rand()/(RAND_MAX + 1.0) * n_nodes);
         do{
            b = (int)(rand()/(RAND_MAX + 1.0) * n_nodes);
         }while(a == b);

         /* calculate current partial energy */
         curr_part_energy = getPartEnergy(matrix, graph_order, a, b, size_column, ones_list, n_nodes, options.alpha);

         /* swap chosen nodes */
         swap(graph_order, a, b);

         /* calculate new partial energy */
         new_part_energy = getPartEnergy(matrix, graph_order, a, b, size_column, ones_list, n_nodes, options.alpha);

         /* calculate new energy */
         new_energy = curr_energy - curr_part_energy + new_part_energy;

         /* generate delta */
         delta_energy = new_energy - curr_energy;

         /* if the new energy is better accept and assign */
         if(delta_energy <= 0){
            curr_energy = new_energy;
            swap_cnt++;
         }else{
            /* if new distance is worse accept but with a probability level */
            expo = expl(-delta_energy/temperature);

            if(expo > (rand()/RAND_MAX)){
               curr_energy = new_energy;
               swap_cnt++;
            }else{
               swap(graph_order, a, b);
            }
         }
      }
      fprintf(energy,"%d\t%.f\t%.3lE\t%.lf\n",step, curr_energy*2, temperature*2, swap_cnt);

      /* cooling process  */
      if(step%options.cooling_interval == 0)
         temperature = temperature * options.cooling_factor;
      progressBar(step,options.mc_steps);
      fflush(stdout);
   }

   printf("\nFINAL ENERGY: %.f\n",curr_energy*2);
   /* save final order  */
   printf("Saving final order...\n");
   fclose(energy);
   if(strrchr(dat_name,'.') == NULL)
      sprintf(order_path,"%s/ordering_%s.dat",dir_name,dat_name);
   else
      sprintf(order_path,"%s/ordering_%s",dir_name,dat_name);
   order = fopen(order_path,"w+");
   if(order == NULL){
      printf("ERROR! File order.dat could not be created.\n");
      exit(1);
   }else{
      fprintf(order,"Protein\tdim1\n");
      for(n=0;n<n_nodes;n++)
         fprintf(order,"%s\t%d\n",nodes[graph_order[n].cur],n);
   }
   fclose(order);

   for(n=0;n<n_nodes;n++){
      free(nodes[n]);
      free(matrix[n]);
      free(ones_list[n]);
   }
   free(nodes);
   free(matrix);
   free(ones_list);
   free(graph_order);

   printf("Done!\n");
   return 0;
}
