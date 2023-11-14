#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

pthread_mutex_t mutexsum; 

typedef struct
{
    int thread_id;
    int start_toss;
    int end_toss;
    long long int *hit_number;
} Arg_toss;

void* calculate(void *args){
    Arg_toss *data = (Arg_toss *)args;

    int thread_id = data->thread_id;
    int start_toss = data->start_toss;
    int end_toss = data->end_toss;
    long long int *hitnumber = data->hit_number;

    long long int local_hit_number = 0;

    unsigned int seed = 8964 + thread_id;

    for (int i = start_toss; i < end_toss ;i++ ){

        double x = ((double) rand_r(&seed) / RAND_MAX) ;
        double y = ((double) rand_r(&seed) / RAND_MAX) ;

        double distance = x*x + y*y;
        
        if (distance < 1.0){
            local_hit_number +=1;
        }
    }

    
    pthread_mutex_lock(&mutexsum);
    *hitnumber += local_hit_number;
    pthread_mutex_unlock(&mutexsum);

    pthread_exit((void *)0);
}

int main(int argc, char *argv[])
{
    int numofthread = atoi(argv[1]);
    long long int numoftoss = atoll(argv[2]);
    int pertoss = (numoftoss/ numofthread );

    pthread_t threads[numofthread];
    


    pthread_mutex_init(&mutexsum, NULL);

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);  

    long long int *hit_number = (long long int *) malloc(sizeof(*hit_number));
    *hit_number = 0; 

    Arg_toss arg[numofthread];
    for (int i =0; i < numofthread ; i++){
        arg[i].thread_id = i;
        arg[i].start_toss = pertoss*i;
        arg[i].end_toss = pertoss*(i+1);
        arg[i].hit_number = hit_number;
        
        pthread_create(&threads[i],&attr,calculate,(void *) &arg[i]);
    }
    pthread_attr_destroy(&attr);

    void *state;
    for (int i = 0;i < numofthread;i++) {
        pthread_join(threads[i], &state);
    }

    pthread_mutex_destroy(&mutexsum);

    double pi = 4 * ((*hit_number) / (double)numoftoss);
    printf("%f\n", pi);

    return 0;    

}
