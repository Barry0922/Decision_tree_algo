#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define MAX_LINE_SIZE 20
#define feature_size 4
#define class_number 3
#define class_name 10
#define sample_number 150
#define random_size 0.7

float feature_importance[feature_size] = {0};
float feature_importance_entropy[feature_size] = {0};
char *feature_name[4] = { "sepal length (cm)", "sepal width (cm)", "petal length (cm)", "petal width (cm)"};
char *classes_name[3] = { "Setosa", "Versicolor", "Virginica"};

typedef struct decision_tree
{
    int node[2];
    int th_level;
    float value;
    float entropy;
    int sample_distribution[class_number];
    int child_node_samples[2];
    struct decision_tree *f_false;
    struct decision_tree *t_true;
}node;

void print_rank(float** arr, int all_num)
{
    for(int i=0; i <all_num ; i++)
    {
        printf(" arr[%d] = %4.2f %4.2f %4.2f %4.2f %4.2f\n", i, arr[i][0], arr[i][1], arr[i][2], arr[i][3], arr[i][4]);
    }
}

void print_rank_2(float** arr, int all_num)
{
    for(int i=0; i <all_num ; i++)
    {
        printf(" arr[%d] = %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n", i, arr[i][0], arr[i][1], arr[i][2], arr[i][3], arr[i][4], arr[i][5]);
    }
}

void shuffle(int *array, int n, int num_shuffles) {
    srand((unsigned)time(NULL));
    for (int j = 0; j < num_shuffles; j++) {
        for (int i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

int* creat_randon(void )
{
    int *random_position;
    random_position = (int*)malloc(sample_number*sizeof(int));
    
    for(int i =0; i < sample_number; i++){
        random_position[i] = i;
    }
    shuffle(random_position, sample_number, sample_number);
    
    /*for(int i =0; i < sample_number; i++){
        printf("%d %d\n", i, random_position[i]);
    }*/
    return random_position;
}

float calculate_entropy(int *p_for_each_classes, int n)
{
    float ans=0.0, p=0.0;
    int count = 0;

    for(int i =0; i <n;i++)
    {
        count = count + *(p_for_each_classes+i);
    }
    for(int i =0; i < n ; i++)
    {
        p = p_for_each_classes[i] / (float)count;
        //printf("%d / %4.2f = %4.2f\n",p_for_each_classes[i], (float)count, p );
        
        if (p==1.0)
        {
            ans = ans +0.0;
        } 
        else if (p!=0.0)
        {
            ans = ans + (p*log2f(p));
        } 
    }
    ans = ans * -1;
    //printf("%4.4f \n", ans);
    return ans;
}

float calculate_gini(int *p_for_each_classes, int n)
{
    float ans=0.0, p=0.0;
    int count = 0;

    for(int i =0; i <n;i++)
    {
        count = count + *(p_for_each_classes+i);
    }
    for(int i =0; i < n ; i++)
    {
        p = p_for_each_classes[i] / (float)count;
        //printf("%d / %4.2f = %4.2f\n",p_for_each_classes[i], (float)count, p );
        
        ans = ans + (p*p);
    }
    ans = 1.0 - ans ;
    //printf("%4.4f \n", ans);
    return ans;
}

void quick_sort_2d_pointers(float **arr, int first_index, int last_index, int th_item, int number_of_items) 
{
    /*
        using quick sor method to sort 2d pointer array (pointer of pointer)

        first_index = index of first *arr
        last_index = index of last *arr
        th_item to sort by which valuse  >> *(arr+ th_item)
        number_of_itmes means how many data (i) in *(arr+i)
    */
    int pivotIndex, index_a, index_b;
    float *temp = malloc(number_of_items*sizeof(float));

    if (first_index < last_index) 
    {
        pivotIndex = first_index;
        index_a = first_index;
        index_b = last_index;

        while (index_a < index_b) 
        {
            while (arr[index_a][th_item] <= arr[pivotIndex][th_item] && index_a < last_index) 
            {
                index_a++;
            }
            
            while (arr[index_b][th_item] > arr[pivotIndex][th_item]) 
            {
                index_b--;
            }

        if (index_a < index_b) 
        {
            temp = *(arr+index_a);
            //printf("%d %d\n", *(temp+0), *(temp+1));
            *(arr+index_a) = *(arr+index_b);
            *(arr+index_b) = temp;
        }
    }

    temp = *(arr+pivotIndex);
    *(arr+pivotIndex) = *(arr+index_b);
    *(arr+index_b) = temp;

    quick_sort_2d_pointers(arr, first_index, index_b - 1, th_item, number_of_items);
    quick_sort_2d_pointers(arr, index_b + 1, last_index, th_item, number_of_items);   
    }
}

void printf_tree(node* t)
{   
    /*
        just print what data have beend learnt in a tree node
    */
    if (t->entropy != 0.0)
    {
        printf("                                                            split node is [%d, %d] in %d th level tree\n", t->node[0], t->node[1], t->th_level);
        printf("                                                            entropy is %4.2f \n", t->entropy);
        printf("                                                            %s is less than %4.2f \n", feature_name[t->node[1]], t->value);
        printf("                                                            current sample distribution %d, %d, %d \n", t->sample_distribution[0], t->sample_distribution[1], t->sample_distribution[2]);
        printf("                                                            after splitting =  [ less : %d , more : %d ] \n", t->child_node_samples[0], t->child_node_samples[1]);
        printf("                                                            current tree is in  %p\n", t);
    }
    else
    {
        int i=0;
        printf("This node is leaf in %d level because ", t->th_level);
        printf("entropy is %4.2f \n", t->entropy);
        printf("current sample distribution %d, %d, %d \n", t->sample_distribution[0], t->sample_distribution[1], t->sample_distribution[2]);
        while (t->sample_distribution[i]==0.0)
        {
            i =i+1;
        }
        printf("all samples here, are belong to %d ( %s ) \n", i, classes_name[i]);
        printf("current tree is in  %p\n", t);
    }
    
}

int* calculate_distribution(float **array, int sample_n, int class_n)
{
    /* 
    arr[4] is put the real data.
    here, to calculate data distribution
    here only fit 3 classs numbers (0, 1, 2), which means only fit the iris dataset.
    If want to implement other data, needs some revised.
    */
    int *dis;
    dis = (int*)malloc(class_n*sizeof(int));
    memset(dis, 0, class_n*sizeof(int));

    for(int i=0; i < sample_n ; i++)
    {   
        if ( array[i][4] == 0.0 )
        {
            dis[0] = dis[0]+1;
        }
        else if ( array[i][4] == 1.0)
        {
            dis[1] = dis[1] +1;
        }
        else if ( array[i][4] == 2.0 )
        {
            dis[2] = dis[2] +1;
        }
    }
    /*
    for(int j =0; j<class_n; j++)
    {
        printf("classes: %d  number: %d \n", j, dis[j] );
    }*/
    return dis;
}

node* decision_tree_recursion(float **t_data, int s_number, int level, int original_smaples_n)
{
    /*
    Input 2D  " training data" array, and using greedy search to find best features to split "training data" array into sub-data.
    Then, used tree stucture to record all relevant information, which is used to identify testing. In other words, record relevant information stand for the learning processing.
    And using ecursion method to spilit sub-data until the entropy of sub data is 0, whihc means growing the decision tree.
    */
    
    //printf(" %d th level, sample number : %d\n", level, s_number);

    int *current_distribution;
    current_distribution = (int*)malloc( class_number*sizeof(int));
    current_distribution = calculate_distribution(t_data, s_number, class_number);
    
    float current_entropy, current_gini;
    current_entropy = calculate_entropy(current_distribution, class_number);
    current_gini = calculate_gini(current_distribution, class_number);
    //printf("current entropy = %4.2f\n", current_entropy);

    if (current_entropy != 0.0)
    {
        float information_gain = -99, target=0.0, tem_true_entropy=0.0, tem_false_entropy=0.0;
        int node_j=0, node_i=0, f_total=0, t_total=0, final_f_total=0, final_t_total=0;

        for(int i =0 ; i < feature_size; i++)
        {
            quick_sort_2d_pointers(t_data, 0, (s_number-1), i, (feature_size+1));

            for(int j =0 ; j < s_number; j++)
            {
                int true_class[3] = {0,0,0}, false_class[3]={0,0,0};
                f_total=0;
                t_total=0;
                target = t_data[j][i];
            
                for(int k=0; k < s_number ; k++)
                {
                    if (t_data[k][i] <= target)
                    {
                        if (t_data[k][4] == 0.0)
                        {
                            false_class[0] = false_class[0]+1;
                        }
                        else if (t_data[k][4] == 1.0)
                        {
                            false_class[1] = false_class[1]+1;
                        }
                        else
                        {
                            false_class[2] = false_class[2]+1;
                        }
                        f_total = f_total+1;
                    }
                    else
                    {
                        if (t_data[k][4] == 0.0)
                        {
                            true_class[0] = true_class[0]+1;
                        }
                            else if (t_data[k][4] == 1.0)
                        {
                            true_class[1] = true_class[1]+1;
                        }
                        else
                        {
                            true_class[2] = true_class[2]+1;
                        }
                        t_total= t_total +1;
                    }
                }
                
                tem_false_entropy = calculate_entropy(false_class, class_number);
                tem_true_entropy = calculate_entropy(true_class, class_number);
                
                if ( ( current_entropy - ( ( tem_false_entropy * f_total/(s_number) ) + ( tem_true_entropy*(t_total)/(s_number) ) ) ) >= information_gain )
                {
                    information_gain = ( current_entropy - ( ( tem_false_entropy*f_total/(s_number) ) + ( tem_true_entropy*(t_total)/(s_number) ) ) );
                    node_j = j; //th sample
                    node_i = i; //th class 
                    final_f_total = f_total;
                    final_t_total = t_total;
                    //printf(" %d, %d ,t: %d %d %d %4.2f ",i,j, true_class[0], true_class[1], true_class[2], tem_true_entropy);
                    //printf("f: %d %d %d %4.2f ", false_class[0], false_class[1], false_class[2], tem_false_entropy);
                    //printf("ig = %4.2f, f_total = %d t_total = %d\n", information_gain, f_total, t_total );
                }
            }
        } 
        
        t_total = final_t_total;
        f_total = final_f_total;

        quick_sort_2d_pointers(t_data, 0, (s_number-1), node_i, feature_size+1);
        //print_rank(t_data, s_number);
        //printf(" node : %d , %d,  feature class : %d feature value : %4.2f t: %d f: %d\n", node_j, node_i, node_i, t_data[node_j][node_i], t_total, f_total);
        //float **ff_group, **tt_group;    
        //ff_group = (float**)malloc((f_total)*sizeof(float*));
        //tt_group = (float**)malloc(((t_total))*sizeof(float*));
    
        float *ff_group[f_total], *tt_group[t_total]; 

        for(int m = 0; m< f_total; m++)
        {
            ff_group[m] = t_data[m];   
        }

        for(int n = f_total ; n < s_number; n++)
        {
            tt_group[n - (f_total)] = t_data[n];
        }

        node *tree, *f_tree, *t_tree;
    
        tree = (node*)malloc(1*sizeof(node));
        f_tree = (node*)malloc(1*sizeof(node));
        t_tree = (node*)malloc(1*sizeof(node));
        
        tree->node[0] = node_j;
        tree->th_level = level;
        tree->node[1] = node_i;
        tree->entropy = current_entropy;
        tree->value = (t_data[node_j][node_i] + t_data[node_j+1][node_i])/2.0;
        tree->child_node_samples[0] = f_total;
        tree->child_node_samples[1] = t_total;
    
        int *right_distribution, *left_distribution;
        float right_entropy, left_entropy, right_gini, left_gini;
    
        right_distribution = (int*)malloc( class_number*sizeof(int));
        left_distribution = (int*)malloc( class_number*sizeof(int));
        right_distribution = calculate_distribution(tt_group, t_total, class_number);
        left_distribution = calculate_distribution(ff_group, f_total, class_number);
    
        right_entropy = calculate_entropy(right_distribution, class_number);
        left_entropy = calculate_entropy(left_distribution, class_number);
        right_gini = calculate_gini(right_distribution, class_number);
        left_gini = calculate_gini(left_distribution, class_number);

        feature_importance[node_i] = feature_importance[node_i]+ ( (current_gini*(float)s_number/(float)original_smaples_n-(float)f_total/(float)original_smaples_n*left_gini - (float)t_total/(float)original_smaples_n*right_gini)/(float)original_smaples_n );
        feature_importance_entropy[node_i] = feature_importance_entropy[node_i] + ( (current_entropy*(float)s_number/(float)original_smaples_n-(float)f_total/(float)original_smaples_n*left_entropy - (float)t_total/(float)original_smaples_n*right_entropy)/(float)original_smaples_n );
        //feature_importance[node_i] = feature_importance[node_i]+ ( (current_entropy*(float)s_number/(float)original_smaples_n-(float)f_total/(float)original_smaples_n*left_entropy - (float)t_total/(float)original_smaples_n*right_entropy)/(float)original_smaples_n );
        //printf("right en :　%4.2f [%d %d %d] , left en %4.2f  [%d %d %d]\n", right_entropy, right_distribution[0], right_distribution[1], right_distribution[2], left_entropy, left_distribution[0], left_distribution[1], left_distribution[2]);
    
        for(int i =0; i < class_number; i++)
        {
            tree->sample_distribution[i] = *(current_distribution+i);
        
        }
        
        free(current_distribution);
        free(left_distribution);
        free(right_distribution);

        printf("                                                            ------------------------ node ------------------------\n");
        printf_tree(tree);
        printf("                                                            ------------------------------------------------------\n\n");

        level = level +1;
    
        f_tree = decision_tree_recursion(ff_group, f_total, level, original_smaples_n);
        t_tree = decision_tree_recursion(tt_group, t_total, level, original_smaples_n);

        tree->f_false = f_tree;
        tree->t_true = t_tree;

        return tree;
    }
    else
    {
        node *tree;
        tree = (node*)malloc(1*sizeof(node));
        tree->node[0] = 0;
        tree->node[1] = 0;
        tree->entropy = current_entropy;
        tree->th_level = level;
        
        for(int i =0; i < class_number; i++)
        {
            tree->sample_distribution[i] = *(current_distribution+i);
        }

        tree->value = 0.0;
        tree->child_node_samples[0]= 0;
        tree->child_node_samples[1] = 0;
        tree->f_false = NULL;
        tree->t_true = NULL;
        
        printf("************************ leaf ************************\n");
        printf_tree(tree);
        printf("******************************************************\n\n");
        
        free(current_distribution);
        
        return tree;
    }
}

float* implement_decision_tree(node* tree, float **t_data, int s_number)
{
    /*
    All rules has been learnt and saved into tree structure.
    By using DFS, follow tree rules to split "testing data" array
    testing data[4] is real answer, testing data[5] is predicted answer
    */
    printf(" tree address %p, data address %p, sample number : %d\n", tree, t_data, s_number);

    int *current_distribution = malloc( class_number*sizeof(int) );
    float current_entropy, target=0.0;
    int node_j=0, node_i=0, f_total=0, t_total=0, ans;

    target = tree->value;
    node_j = tree->node[0];
    node_i = tree->node[1];
    current_entropy = tree->entropy;

    if (current_entropy != 0.0)
    {  
        quick_sort_2d_pointers(t_data, 0, (s_number-1), node_i, feature_size+2);
        
        for(int k=0; k < s_number ; k++)
        {
            if (t_data[k][node_i] <= target)
            {
                f_total = f_total+1;
            }
            else
            {    
                t_total= t_total +1;
            }
            
        }
        printf(" node : %d , %d \n", node_j, node_i);
        printf(" feature class : %4.2f  ( %s ) \n", target ,feature_name[node_i]);
        printf(" after splitting, [ f : %d, t : %d] \n\n", f_total, t_total);
    
        float *ff_group[f_total], *tt_group[t_total];
        for(int m = 0; m< f_total; m++)
        {
            ff_group[m] = t_data[m];
        }
        
        for(int n = f_total ; n < s_number; n++)
        {
            tt_group[n - (f_total)] = t_data[n];
        }

        if (tree->f_false != NULL)
        {
            implement_decision_tree(tree->f_false, ff_group, f_total);
        }

        if (tree->t_true != NULL)
        {
            implement_decision_tree(tree->t_true, tt_group, t_total);    
        }
        free(tree);
    }
    else
    {
        current_distribution = tree->sample_distribution;
        for(int i=0; i < class_number; i++)
        {
            if (current_distribution[i] != 0.0)
            {
                ans = i;
                break;
            }
        }
        
        printf(" Because entropy is 0, this sample belongs to %s (%d) \n", classes_name[ ans ], ans);
        for(int i=0; i < s_number; i++)
        {
            t_data[i][class_number+2] = (float)ans;
            printf(" Sample ( %4.2f, %4.2f, %4.2f, %4.2f  real ans : %s (%4.2f) ) in addr %p belongs to %s (%4.2f)\n", t_data[i][0], t_data[i][1], t_data[i][2], t_data[i][3], classes_name[ (int)t_data[i][4] ], t_data[i][4], t_data[i], classes_name[(int)t_data[i][5]], t_data[i][5] );            
        }
        printf("\n");
        free(tree);
    }
}

float calculate_score(float **arr, int s_number, int samples, int ans)
{
    /*
    since arr[4] is real data, ans[5] is ml predicted data
    compare above data to calculate the final accuracy 
    */
    float score, right=0.0, wrong=0.0;

    for(int i =0; i <s_number ; i++)
    {   
        //printf("%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f in %p \n", arr[i][0], arr[i][1], arr[i][2], arr[i][3], arr[i][samples], arr[i][ans], arr[i] );
        if (arr[i][samples] == arr[i][ans])
        {
            right = right+1;
        } 
        else
        {
            wrong = wrong +1;
        }
    }
    score = right/(wrong+right);
    return score;
}

void print_feature_importantce(float *arr)
{
    float total =0;
    for(int i =0 ; i < feature_size; i++)
    {
        total = total + arr[i];
    }
    printf("Feature importances : \n");
    printf(" [ ");
    for(int i =0 ; i < feature_size; i++)
    {
        printf("%4.4f ", arr[i]/total);
    }
    printf(" ]\n");
}

int main(void) 
{
    //---------------------------------- this part is used to featch the iris data into "rank array" from CSV file------------------------------------------------------------------------------
    char *file_name = "iris_03.csv";
    FILE *fp;
    fp = fopen(file_name, "r");

    if (!fp) 
    {
        fprintf(stderr, "failed to open file for reading\n");
    }
    
    char line[MAX_LINE_SIZE], *result = NULL;
    
    float **rank;
    rank = (float**)malloc(sample_number* sizeof(float*));
    for(int i = 0; i <sample_number ; i++)
    {
        *(rank+i) = (float*)malloc((feature_size+1) * sizeof(float));
    }

    int th_sample = 0;
    while(fgets(line, MAX_LINE_SIZE, fp) != NULL)
    {
        
        result = strtok(line, ",");
        int i=0;

        while( result != NULL ) 
        {       
            rank[th_sample][i]= atof(result);
            result = strtok(NULL, ",");
            i++;
        }   
        th_sample = th_sample +1;
    }
    //print_rank(rank, sample_number);
    //quick_sort_2d_pointers(rank, 0, (sample_number-1), 0, (feature_size+1) );
    //print_rank(rank, sample_number);

    fclose (fp);
    
    //---------------------------------- split data from rank array into training and testing data array ------------------------------------------------------------------------------
    int *random, *temp; 
    float **training_data, **testing_data;
    
    testing_data = (float**)malloc((sample_number*(1-random_size))* sizeof(float*));
    training_data = (float**)malloc((sample_number*random_size)* sizeof(float*));

    for(int i = 0; i <sample_number*random_size ; i++)
    {
        training_data[i] = (float*)malloc((feature_size+1) * sizeof(float));
    }

    for(int i = 0; i <(sample_number*(1-random_size)) ; i++)
    {
        testing_data[i] = (float*)malloc((feature_size+2) * sizeof(float));
    }
    /*// -----------------------------------------------for case where doesn't use fixed random state-----------------------------------------------------------------------------------------
    random = creat_randon();
    for(int i =0; i < sample_number*random_size; i++)
    {
        training_data[i][0] = rank[random[i]][0];
        training_data[i][1] = rank[random[i]][1];
        training_data[i][2] = rank[random[i]][2];
        training_data[i][3] = rank[random[i]][3];
        training_data[i][4] = rank[random[i]][4];
    }

    int a = (int)(sample_number*random_size);
    for(int i = a ; i < sample_number; i++)
    {
        testing_data[i-a][0] = rank[random[i]][0];
        testing_data[i-a][1] = rank[random[i]][1];
        testing_data[i-a][2] = rank[random[i]][2];
        testing_data[i-a][3] = rank[random[i]][3];
        testing_data[i-a][4] = rank[random[i]][4];
        testing_data[i-a][5] = 9.9;    
    }
    free(random);
    */

    // --------------------------------- split data in training and testing by ratio (random size)---------------------------------
    // --------------------------------- testing[i][5] 5th columns is predicted answer---------------------------------------------
    for(int i =0; i < sample_number*random_size; i++)
    {
        training_data[i][0] = rank[i][0];
        training_data[i][1] = rank[i][1];
        training_data[i][2] = rank[i][2];
        training_data[i][3] = rank[i][3];
        training_data[i][4] = rank[i][4];
    }

    int a = (int)(sample_number*random_size);
    for(int i = a ; i < sample_number; i++)
    {
        testing_data[i-a][0] = rank[i][0];
        testing_data[i-a][1] = rank[i][1];
        testing_data[i-a][2] = rank[i][2];
        testing_data[i-a][3] = rank[i][3];
        testing_data[i-a][4] = rank[i][4];
        testing_data[i-a][5] = 9.9;    
    }

    //print_rank(training_data, (sample_number*random_size));
    quick_sort_2d_pointers(training_data, 0, ((sample_number*random_size)-1), 0, (feature_size+1) );    
    //print_rank(training_data, (sample_number*random_size));
    //printf("---------------above is training -----------below is testing -----------\n");
    //print_rank_2(testing_data, (sample_number*(1-random_size)) );
    quick_sort_2d_pointers(testing_data, 0, ((sample_number*(1-random_size))-1), 0, (feature_size+2) );    
    //print_rank_2(testing_data, (sample_number*(1-random_size)) );
    
    //------------------------ create decistion node stucture to save node data ----------------------------------------------------------
    node *decision_tree_nodes;
    decision_tree_nodes = (node*)malloc(sizeof(node));
    
    decision_tree_nodes =  decision_tree_recursion( training_data, (sample_number*random_size), 0, (sample_number*random_size) );
    printf("\n--------------------------------above are learning process--------------------------------\n--------------------------------below are predicting process--------------------------------\n\n");
    implement_decision_tree(decision_tree_nodes, testing_data, (sample_number*(1-random_size)));
    
    float final_accuracy;
    final_accuracy = calculate_score(testing_data, (sample_number*(1-random_size)), 4, 5 );

    printf("final accuracy = %4.4f \n", final_accuracy);
    printf("\n------------------------------\n");
    print_feature_importantce(feature_importance_entropy);
    printf("--------------------------------\n");
    //print_feature_importantce(feature_importance);

    free(rank);
    free(training_data);
    free(testing_data);
    for(int i =0; i < sample_number*random_size; i++)
    {
        free( *(training_data+i) );
        free( *(rank+i) );
    }
    int b = (int)(sample_number*random_size);
    for(int i = b ; i < sample_number; i++)
    {
        free( *(testing_data+(i-a)) ); 
        free( *(rank+i) );    
    }
    free(decision_tree_nodes);

    return 0;
}