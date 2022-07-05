/* 
 * G01
 * Integrantes:
 * - Jhon Steeven Cabanilla Alvarado
 * - Miguel Chaveinte García
 * 
 * Probabilistic approach to locate maximum heights
 * Hill Climbing + Montecarlo
 *
 * MPI version
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2021/2022
 *
 * v1.1
 *
 * (c) 2022 Arturo Gonzalez Escribano
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<limits.h>
#include<sys/time.h>
#include<omp.h>

/* Headers for the MPI assignment versions */
#include<mpi.h>                              
#include<stddef.h>                           

/* 
 * Global variables and macro to check errors in calls to MPI functions
 * The macro shows the provided message and the MPI string in case of error
 */
char mpi_error_string[ MPI_MAX_ERROR_STRING ];
int mpi_string_len;
#define MPI_CHECK( msg, mpi_call )	{ int check = mpi_call; if ( check != MPI_SUCCESS ) { MPI_Error_string( check, mpi_error_string, &mpi_string_len); fprintf(stderr,"MPI Error - %s - %s\n", msg, mpi_error_string ); MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE ); } }


#define	PRECISION	10000
#define min(X, Y) (((X) < (Y)) ? (X) : (Y))  //MODIFICACION

/* 
 * Structure to represent a climbing searcher 
 * 	This structure can be changed and/or optimized by the students
 */
typedef struct {
	int id;			// Searcher identifier
	int pos_row, pos_col;		// Position in the grid
	int steps;			// Steps count
	int follows;			// When it finds an explored trail, who searched that trail
} Searcher;


/* 
 * Function to get wall time
 */
double cp_Wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}

/* 
 * Macro function to simplify accessing with two coordinates to a flattened array
 * 	This macro-function can be changed and/or optimized by the students
 */
#define accessMat( arr, exp1, exp2 )	arr[ (int)(exp1) * columns + (int)(exp2) ]

#define get_height( x, y, rows, columns, x_min, x_max, y_min, y_max ) (int)( PRECISION * (float)(2 * sin((float)(x_min) + ( ((float)(x_max) - (float)(x_min)) / (int)(rows) ) * (int)(x)) * cos(((float)(y_min) + ( ((float)(y_max) - (float)(y_min)) / (int)(columns) ) * (int)(y))/2) + log( fabs(((float)(y_min) + ( ((float)(y_max) - (float)(y_min)) / (int)(columns) ) * (int)(y)) - M_PI_2) )) );

/*
 * Function: Generate height for a given position
 * 	This function can be changed and/or optimized by the students

int get_height( int x, int y, int rows, int columns, float x_min, float x_max, float y_min, float y_max  ) {
	//Calculate the coordinates of the point in the ranges 
	float x_coord = x_min + ( (x_max - x_min) / rows ) * x;
	float y_coord = y_min + ( (y_max - y_min) / columns ) * y;
	//Compute function value 
	float value = 2 * sin(x_coord) * cos(y_coord/2) + log( fabs(y_coord - M_PI_2) );
	//Transform to fixed point precision 
	int fixed_point = (int)( PRECISION * value );
	return fixed_point;
} */

/*
 * Function: Climbing step
 * 	This function can be changed and/or optimized by the students
 */
int climbing_step( int rows, int columns, Searcher *searchers, int search, int *heights, int *trails, float x_min, float x_max, float y_min, float y_max, int ini, int fin) {
	int search_flag = 0;

	/* Annotate one step more, landing counts as the first step */
	searchers[ search ].steps ++;

	/* Get starting position */
	int pos_row = searchers[ search ].pos_row; 
	int pos_col = searchers[ search ].pos_col;
	
	int pos_row_local = pos_row-ini; //MODIFICACION: restamos para obtener la coordenada local

	/* Stop if searcher finds another trail */
	int check;

	check = accessMat( trails, pos_row_local, pos_col );
	accessMat( trails, pos_row_local, pos_col ) = searchers[ search ].id;
	
	


	if ( check != -1 ) {
		accessMat( trails, pos_row_local, pos_col ) = check;
		search_flag = 1;
		searchers[ search ].follows = check;
	}
	else {
		/* Annotate the trail */
		//accessMat( trails, pos_row_local, pos_col ) = searchers[ search ].id;

		/* Locate the highest climbing direction */
		float local_max = accessMat( heights, pos_row_local, pos_col );
		int climbing_direction = 0;
		float altura=0;
		
		if ( pos_row > 0 ) {
			
			if (pos_row_local==0){
				altura = accessMat( heights, fin-ini, pos_col); 
			 	if(altura==INT_MIN) altura = accessMat( heights, fin-ini, pos_col )=get_height( pos_row-1, pos_col, rows, columns, x_min, x_max, y_min, y_max );
			}
			else{
				altura = accessMat( heights, pos_row_local-1, pos_col );
				if(altura==INT_MIN) altura = accessMat( heights, pos_row_local-1, pos_col )=get_height( pos_row-1, pos_col, rows, columns, x_min, x_max, y_min, y_max );
			}
			 
			
			/* Annotate the travelling direction if higher */
			if ( altura > local_max ) {
				climbing_direction = 1;
				local_max = altura;
			}
		}
		
		
		if ( pos_row < rows-1 ) {
			
			
			if (pos_row_local==fin-ini-1){
				altura = accessMat( heights, fin-ini+1, pos_col); 
			 	if(altura==INT_MIN) altura = accessMat( heights, fin-ini+1, pos_col )=get_height( pos_row+1, pos_col, rows, columns, x_min, x_max, y_min, y_max );
			}
			else{
				altura = accessMat( heights, pos_row_local+1, pos_col );
				if(altura==INT_MIN) altura = accessMat( heights, pos_row_local+1, pos_col )=get_height( pos_row+1, pos_col, rows, columns, x_min, x_max, y_min, y_max );
			}
			
			
			
			/* Annotate the travelling direction if higher */
			if ( altura > local_max ) {
				climbing_direction = 2;
				local_max = altura;
			}
		}
		
		
		if ( pos_col > 0 ) {
			altura = accessMat( heights, pos_row_local, pos_col-1 );
			
			if(altura==INT_MIN) altura = accessMat( heights, pos_row_local, pos_col-1 )=get_height( pos_row, pos_col-1, rows, columns, x_min, x_max, y_min, y_max );
							
			/* Annotate the travelling direction if higher */
			if ( altura > local_max ) {
				climbing_direction = 3;
				local_max = altura;
			}
		}
		
		
		if ( pos_col < columns-1 ) {
			altura = accessMat( heights, pos_row_local, pos_col+1 );
			if(altura==INT_MIN) altura = accessMat( heights, pos_row_local, pos_col+1 )=get_height( pos_row, pos_col+1, rows, columns, x_min, x_max, y_min, y_max );
			
			/* Annotate the travelling direction if higher */
			if ( altura > local_max ) {
				climbing_direction = 4;
				local_max = altura;
			}
		}
		
		
		
		if(climbing_direction==0){
			searchers[ search ].follows = searchers[ search ].id;
			search_flag = 1;
		}
		else if(climbing_direction==1) pos_row--;
		else if (climbing_direction==2) pos_row++;
		else if (climbing_direction==3) pos_col--;
		else pos_col++;
		
		searchers[ search ].pos_row = pos_row; 
		searchers[ search ].pos_col = pos_col;
		
		/*MODIFICACION*/
		if(ini>pos_row ){
			search_flag = 2;
			//printf("por abajo");
		}
		
		if(pos_row >=fin){
			search_flag = 3;
			//printf("por arriba");
		}
	}

	/* Return a flag to indicate if search should stop */
	return search_flag;
}


#ifdef DEBUG
/* 
 * Function: Print the current state of the simulation 
 */
void print_heights( int rows, int columns, int *heights ) {
	/* 
	 * You don't need to optimize this function, it is only for pretty 
	 * printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i,j;
	printf("Heights:\n");
	printf("+");
	for( j=0; j<columns; j++ ) printf("-------");
	printf("+\n");
	for( i=0; i<rows; i++ ) {
		printf("|");
		for( j=0; j<columns; j++ ) {
			char symbol;
			if ( accessMat( heights, i, j ) != INT_MIN ) 
				printf(" %6d", accessMat( heights, i, j ) );
			else
				printf("       ");
		}
		printf("|\n");
	}
	printf("+");
	for( j=0; j<columns; j++ ) printf("-------");
	printf("+\n\n");
}

void print_trails( int rows, int columns, int *trails ) {
	/* 
	 * You don't need to optimize this function, it is only for pretty 
	 * printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i,j;
	printf("Trails:\n");
	printf("+");
	for( j=0; j<columns; j++ ) printf("-------");
	printf("+\n");
	for( i=0; i<rows; i++ ) {
		printf("|");
		for( j=0; j<columns; j++ ) {
			char symbol;
			if ( accessMat( trails, i, j ) != -1 ) 
				printf("%7d", accessMat( trails, i, j ) );
			else
				printf("       ", accessMat( trails, i, j ) );
		}
		printf("|\n");
	}
	printf("+");
	for( j=0; j<columns; j++ ) printf("-------");
	printf("+\n\n");
}
#endif // DEBUG

/*
 * Function: Print usage line in stderr
 */
void show_usage( char *program_name ) {
	fprintf(stderr,"Usage: %s ", program_name );
	fprintf(stderr,"<rows> <columns> <x_min> <x_max> <y_min> <y_max> <searchers_density> <short_rnd1> <short_rnd2> <short_rnd3>\n");
	fprintf(stderr,"\n");
}


/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
	// This eliminates the buffer of stdout, forcing the messages to be printed immediately
	setbuf(stdout,NULL);

	int i,j;

	// Simulation data
	int rows, columns;		// Matrix sizes
	float x_min, x_max;		// Limits of the terrain x coordinates
	float y_min, y_max;		// Limits of the terrain y coordinates

	float searchers_density;	// Density of hill climbing searchers
	unsigned short random_seq[3];	// Status of the random sequence

	int *heights = NULL;		// Heights of the terrain points
	int *trails = NULL;		// Searchers trace and trails
	//int *tainted = NULL;		// Position found in a search
	int *follows = NULL;		// Compacted list of searchers "follows"
	int num_searchers;		// Number of searchers
	int num_searchers_total;	//MODIFICACION
	Searcher *searchers = NULL;	// Searchers data
	int *total_steps = NULL;	// Annotate accumulated steps to local maximums

	/* 0. Initialize MPI */
	MPI_Init( &argc, &argv );
	int rank;
	int nprocs; /*MODIFICACION*/
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &nprocs ); /*MODIFICACION*/
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

	/* 1. Read simulation arguments */
	/* 1.1. Check minimum number of arguments */
	if (argc != 11) {
		fprintf(stderr, "-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 1.2. Read argument values */
	rows = atoi( argv[1] );
	columns = atoi( argv[2] );
	x_min = atof( argv[3] );
	x_max = atof( argv[4] );
	y_min = atof( argv[5] );
	y_max = atof( argv[6] );
	searchers_density = atof( argv[7] );

	/* 1.3. Read random sequences initializer */
	for( i=0; i<3; i++ ) {
		random_seq[i] = (unsigned short)atoi( argv[8+i] );
	}


#ifdef DEBUG
	/* 1.4. Print arguments */
	if ( rank == 0 ) {
		printf("Arguments, Rows: %d, Columns: %d\n", rows, columns);
		printf("Arguments, x_range: ( %d, %d ), y_range( %d, %d )\n", x_min, x_max, y_min, y_max );
		printf("Arguments, searchers_density: %f\n", searchers_density );
		printf("Arguments, Init Random Sequence: %hu,%hu,%hu\n", random_seq[0], random_seq[1], random_seq[2]);
		printf("\n");
	}
#endif // DEBUG


	/* 2. Start global timer */
	MPI_CHECK( "Clock: Start-Barrier ", MPI_Barrier( MPI_COMM_WORLD ) );
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */	
 
	///*Implementacion Struct Searcher*///
	int fields = 5;
	
	int array_of_blocklengths[] = {1,1,1,1,1};
	
	// Block displacements
	MPI_Aint array_of_displacements[] = {
		offsetof( Searcher, id ),
		offsetof( Searcher, pos_row ),
		offsetof( Searcher, pos_col),
		offsetof( Searcher, steps ),
		offsetof( Searcher, follows) 
		};
		
	// Block types
	MPI_Datatype array_of_types[] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
	
	MPI_Datatype MPI_Searcher;
	
	MPI_CHECK( "creacion MPI_Searcher",MPI_Type_create_struct(fields, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_Searcher ));

	MPI_CHECK( "subida MPI_Searcher",MPI_Type_commit( &MPI_Searcher));
	
		
	
 	//Comenzamos realizando la division por procesos
 	int resto = rows % nprocs;
 	int my_size = rows / nprocs; 	//Tamaño de cada uno de los procesos
 	int my_begin = rank * my_size;	//Comienzo de cada uno de los procesos
 	
 	
 	
 	//Caso en el que el size no sea divisible entre el numero de procesos
 	if(resto!=0){
	 	if(rank < resto){
	 		my_size += 1;
	 	}
	 	my_begin += min(rank, resto);
 	}
 	
 	
 	
 	
 	
	/* Statistical data */
	int num_local_max = 0;
	int max_accum_steps = INT_MIN;
	int total_tainted = 0;
	int max_height = INT_MIN;
	int maximo_local=INT_MIN; /*MODIFICACION*/
	int recorridos_local=0; /*MODIFICACION*/
	


                         
	/* 3. Initialization */
	/* 3.1. Memory allocation */
	num_searchers = (int)( rows * columns * searchers_density );
	num_searchers_total = num_searchers;

	searchers = (Searcher *)malloc( sizeof(Searcher) * num_searchers ); 
	total_steps = (int *)malloc( sizeof(int) * num_searchers ); 
	follows = (int *)malloc( sizeof(int) * num_searchers );
	if ( searchers == NULL || total_steps == NULL || follows == NULL ) {
		fprintf(stderr,"-- Error allocating searchers structures for size: %d\n", num_searchers );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	
	
	
	//Reservamos solo la parte necesaria de cada matriz
	heights = (int *)malloc( sizeof(int) * (size_t)(my_size+2) * (size_t)columns );
	trails = (int *)malloc( sizeof(int) * (size_t)my_size * (size_t)columns );
	//tainted = (int *)calloc(my_size *columns,sizeof(int));
	

	if ( heights == NULL || trails == NULL ) {
		fprintf(stderr,"-- Error allocating terrain structures for size: %d x %d \n", rows, columns );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	
	

	/* 3.2. Terrain initialization */
	
	//MODIFICACION | Forma de reducir el bucle de INICIALIZACION DEL TERRENO
	for(i=0; i<my_size*columns; i++){
		heights[i] = INT_MIN;
		trails[i] = -1;
		//tainted[i] = 0;
	}
	/*Modificacion añadir halos y rellenar*/
	for(i=my_size; i<(my_size+2)*columns; i++){
		heights[i] = INT_MIN;
	}
	

	
	
	/*3.3 Searchers initialization */
	int search;
	int cont = 0;
	for( search = 0; search < num_searchers; search++ ) {
	
		searchers[ search ].id = search ;
		searchers[ search ].pos_row = (int)( rows * erand48( random_seq ) ); 	//Posicion global
		searchers[ search ].pos_col = (int)( columns * erand48( random_seq ) ); //Posicion global
		searchers[ search ].steps = 0;
		searchers[ search ].follows = -1;
		total_steps[ search ] = 0;
	
		
		if(searchers[ search ].pos_row >= my_begin && searchers[ search ].pos_row < (my_begin + my_size)){ 
		//El searcher pertenece al proceso
			searchers[cont] = searchers[ search ];
			cont ++;
		}
	}
	
	num_searchers = cont;
	
	
	
	//Declaracion arrays -- Comunicaciones
	Searcher *terminados = NULL;		// Searcher terminados
	Searcher *EnviadosDelante = NULL;	// Searchers que no pertenecen
	Searcher *EnviadosDetras = NULL;		// Searchers que no pertenecen

	
	//Memory allocation
	terminados = (Searcher *)malloc( sizeof(Searcher) * num_searchers_total ); 
	EnviadosDelante = (Searcher *)malloc( sizeof(Searcher) * num_searchers_total ); 
	EnviadosDetras = (Searcher *)malloc( sizeof(Searcher) * num_searchers_total );
	
	/*Comprobar que no sean nulos*/
	if ( terminados == NULL || EnviadosDelante == NULL || EnviadosDetras == NULL) {
		fprintf(stderr,"-- Error allocating terrain structures for size: %d x %d \n", rows, columns );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	
	int index_delante = 0;
	int index_detras = 0;
	int index_terminado = 0;
	
		
	/* 4. Compute searchers climbing trails */
	
	int flag_terminados = 1; 
	while(flag_terminados != 0){
	
		for( search = 0; search < num_searchers; search++ ) {
			int search_flag = 0;
			int pos_row=searchers[ search ].pos_row;
			int pos_col = searchers[ search ].pos_col;
			/* Compute the height */
			int pos_row_local = pos_row-my_begin; //MODIFICACION: restamos para obtener la coordenada local
			float altura=accessMat( heights, pos_row_local, pos_col );
			if(altura==INT_MIN) accessMat( heights, pos_row_local, pos_col ) = get_height( pos_row, pos_col, rows, columns, x_min, x_max, y_min, y_max );
			while( ! search_flag ) {
				search_flag = climbing_step( rows, columns, searchers, search, heights, trails, x_min, x_max, y_min, y_max, my_begin, (my_begin + my_size));
			}
			
			//Comprobar posicion de columnas
			if(search_flag==2){
				//Se pasa por debajo
				EnviadosDetras[index_detras] = searchers[ search ];
				//printf("detras");
				index_detras++;
				
			} else if(search_flag==3 ){
				//Se pasa por encima
				EnviadosDelante[index_delante] = searchers[ search ];
				//printf("delante");
				index_delante++;		
			} else{
				//Pertenece
				terminados[index_terminado] = searchers[ search ];
				//printf("terminado");
				index_terminado++;
			}
			
			
		} //Fin For correcto. Obtenemos todos los terminados y procedemos a realizar la comunicacion
			
			
		MPI_Request request1 = MPI_REQUEST_NULL;
		MPI_Request request2 = MPI_REQUEST_NULL;
		
		MPI_Status status1;
		MPI_Status status2;
		
		int searchers_llegados=0;
		int total_recibidos=0;
		
		if(rank != nprocs - 1){
		
			MPI_CHECK("Envios delante menos del ultimo proceso", MPI_Isend(EnviadosDelante, index_delante, MPI_Searcher, rank+1, 1, MPI_COMM_WORLD, &request1) ); /*CHECKS*/
			
			
			MPI_CHECK("Recibo delante menos del utlimo proceso", MPI_Recv(searchers, num_searchers_total, MPI_Searcher, rank+1, 0, MPI_COMM_WORLD, &status1); ); /*CHECKS*/
			
			MPI_CHECK("Obtenemos cantidad de buscadores, proceso 0", MPI_Get_count(&status1, MPI_Searcher, &searchers_llegados) ); /*CHECKS*/
			total_recibidos=total_recibidos+searchers_llegados;
		}
		
		MPI_Wait(&request1, MPI_STATUS_IGNORE);
		index_delante=0;
	
			
		if(rank != 0){
			MPI_CHECK("Envios detras menos del primer proceso", MPI_Isend(EnviadosDetras, index_detras, MPI_Searcher, rank-1, 0, MPI_COMM_WORLD, &request2); ); /*CHECKS*/
			
			
			
			MPI_CHECK("Recibo detras menos del primer proceso", MPI_Recv(&searchers[total_recibidos], num_searchers_total-total_recibidos, MPI_Searcher, rank-1, 1, MPI_COMM_WORLD, &status2) ); /*CHECKS*/
			
			MPI_CHECK("Obtenemos cantidad de buscadores, proceso intermedio", MPI_Get_count(&status2, MPI_Searcher, &searchers_llegados) ); /*CHECKS*/
			total_recibidos=total_recibidos+searchers_llegados;
		}
		
		MPI_Wait(&request2, MPI_STATUS_IGNORE);
		index_detras=0;
		
		
		num_searchers=total_recibidos;
		
		int reduction = 0;
		MPI_CHECK("Reduccion terminados", MPI_Allreduce(&index_terminado, &reduction, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ); /*CHECKS*/
		//printf("iteracion\n");
		
		if(num_searchers_total == reduction){
			flag_terminados = 0;		
		}
				
		
	} //Fin while


#ifdef DEBUG
/* Print computed heights at the end of the search.
 * You can ignore this debug feature.
 * If you want to use this functionality in parallel processes BEWARE: You should 
 * modify it to sincronize the processes to print in order, and each process should
 * prints only its part */
print_heights( rows, columns, heights );
#endif

	/*MODIFICACION*/
	
	 
	
	if(rank != 0){
		MPI_CHECK( "Enviar cont_acabados a proceso 0", MPI_Send(&index_terminado, 1, MPI_INT, 0, 10, MPI_COMM_WORLD) );
		MPI_CHECK( "Enviar acabados a proceso 0", MPI_Send(terminados, index_terminado, MPI_Searcher, 0, 11, MPI_COMM_WORLD) );
	}

	if(rank == 0){
	
		//Searcher *llegados = (Searcher *)malloc( sizeof(Searcher) * num_searchers_total );
	
		int vamos_teniendo=index_terminado;
		//printf("index_termin: %d\n",index_terminado);
				
		for(i = 0; i<vamos_teniendo; i++){
			searchers[i] = terminados[i];
		}
	
		for(i = 1; i<nprocs; i++){
			MPI_Status status1;
			MPI_Status status2;
			
			int llegan=0;
			
			MPI_CHECK( "Recibir index_terminado de cada uno", MPI_Recv(&llegan, 1, MPI_INT, i, 10, MPI_COMM_WORLD, &status1) );
			MPI_CHECK( "Recibir terminados de cada uno", MPI_Recv(&searchers[vamos_teniendo], llegan, MPI_Searcher, i, 11, MPI_COMM_WORLD, &status2) );
			/*for(j = 0; j<llegan; j++){
				searchers[j + vamos_teniendo] = llegados[j];
			}*/
			//printf("llegan: %d\n",llegan);
			vamos_teniendo = vamos_teniendo + llegan;
			
		}
		
		//printf("searche total: %d\n",num_searchers_total);
		//printf("vamos teniendo: %d\n",vamos_teniendo );
		
		
		/* 5. Compute the leading follower of each searcher */
		for( search = 0; search < num_searchers_total; search++ ) {
			follows[ searchers[ search ].id ] = searchers[ search ].follows;
		}
		
		for( search = 0; search < num_searchers_total; search++ ) {
		
			int search_flag = 0;
			int parent = searchers[ search ].id;
			int follows_to = follows[ parent ];
			while( ! search_flag ) {
				if ( follows_to == parent ) search_flag = 1;
				else {
					parent = follows_to;
					follows_to = follows[ parent ];
				}
			}
			searchers[ search ].follows = follows_to;
		}
		
		
		/* 6. Compute accumulated trail steps to each maximum */
		for( search = 0; search < num_searchers_total; search++ ) {
			int pos_max = searchers[ search ].follows;
			total_steps[ pos_max ] += searchers[ search ].steps;
		}
	
		/* 7. Compute statistical data */
		for( search = 0; search < num_searchers_total; search++ ) {
			int id = searchers[ search ].id;
			/* Maximum of accumulated trail steps to a local maximum */
			if ( max_accum_steps < total_steps[ id ] ) 
				max_accum_steps = total_steps[ id ];
			if ( searchers[ search ].follows ==  id ) {
				num_local_max++;
			}
			
		}
		
		//free(llegados);		
	} //FIN COMPUTACION QUE SE NECESITAN CONOCER RANK=0
	
	/* If this searcher found a maximum, check the maximum value */
	for( search = 0; search < index_terminado; search++ ) {
		int id = terminados[ search ].id; // teniamos --> searcher[ search ].id;
		if ( terminados[ search ].follows == id ) {
			int pos_row = terminados[ search ].pos_row - my_begin;
			int pos_col = terminados[ search ].pos_col;
			int valor_altura=accessMat( heights, pos_row, pos_col );
			if ( maximo_local < valor_altura )
				maximo_local = valor_altura;
		}
	}
	//printf("maximo  %d\n",maximo_local);

	MPI_CHECK( "reduce de las alturas maximas", MPI_Reduce(&maximo_local, &max_height, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD) );
	
	for( i=0; i<my_size*columns; i++ ) {
		if ( trails[i] !=-1 ) 
			recorridos_local++;
	}
		
	
	MPI_CHECK( "reduce de casillas manchadas", MPI_Reduce(&recorridos_local, &total_tainted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) );
	
	
/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

	/* 5. Stop global time */
	MPI_CHECK( "End-Barrier", MPI_Barrier( MPI_COMM_WORLD ) );
	ttotal = cp_Wtime() - ttotal;

	/* 6. Output for leaderboard */
	if ( rank == 0 ) { 
		printf("\n");
		/* 6.1. Total computation time */
		printf("Time: %lf\n", ttotal );

		/* 6.2. Results: Statistics */
		printf("Result: %d, %d, %d, %d\n\n", 
			num_local_max,
			max_height,
			max_accum_steps,
			total_tainted );
	}
			
	/* 7. Free resources */	
	free( searchers );
	free( total_steps );
	free( follows );
	free( heights );
	free( trails );
	//free( tainted );
	
	//MODIFICACION
	free( terminados );
	free( EnviadosDelante );
	free( EnviadosDetras );
	MPI_Type_free( &MPI_Searcher);
	

	/* 8. End */
	MPI_Finalize();
	return 0;
}
