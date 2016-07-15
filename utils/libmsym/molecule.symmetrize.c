//
//  msym_example.c
//  libmsym
//
//  Created by Marcus Johansson on 10/09/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "msym.h"

int read_xyz(const char *name, msym_element_t **ratoms) {
    FILE *fp = fopen(name,"r");
    msym_element_t *a;
    int l;
    char buf[1024];
    if(NULL == fp){
        fprintf(stderr, "could not open file: %s\n",name);
        return -1;
    }
    if (NULL == fgets(buf, sizeof(buf), fp) || sscanf(buf," %d ",&l) != 1){
        fprintf(stderr,"Unable to read file %s\n",name);
        fclose(fp);
        return -1;
    }
    if(l < 300000) {
        a = malloc(l*sizeof(msym_element_t));
        memset(a,0,l*sizeof(msym_element_t));
    } else {
        fprintf(stderr, "Too many elements in file %d\n",l);
        fclose(fp);
        return -1;
    }
    
    //char * fgets ( comment, sizeof(comment), fp );
    if(NULL != fgets(buf, sizeof(buf), fp)){
//         printf("Comment: %.*s", (int)sizeof(buf), buf);
    }
    
    for (int i = 0; i < l && fgets(buf, sizeof(buf), fp) && sscanf(buf, "%s %lf %lf %lf", a[i].name, &(a[i].v[0]),  &(a[i].v[1]),  &(a[i].v[2])) == 4 && i < l; i++) {}
    *ratoms = a;
    fclose(fp);
    return l;
    
}

void printSALC(msym_salc_t *salc, msym_element_t *melements){

    
    double (*space)[salc->fl] = (double (*)[salc->fl]) salc->pf;
    for(int d = 0;d < salc->d;d++){
        if(salc->d > 1) printf("Component %d:\n",d+1);
        for(int line = 0; line < salc->fl; line+=6){
            for(int i = line;i < line + 6 && i < salc->fl;i++){
                msym_basis_function_t *bf = salc->f[i];
                printf(" %d%s %-8s\t",(int)(bf->element-melements)+1, bf->element->name,bf->name);
            }
            printf("\n");
            
            for(int i = line;i < line + 6 && i < salc->fl;i++){
                printf("%10.7lf\t", space[d][i]);
                
            }
            printf("\n\n");
        }
        printf("\n");
    }

}

int example(const char* in_file, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *elements = NULL;
    
    const char *error = NULL;
    char point_group[6];
    double cm[3], radius = 0.0, symerr = 0.0;
    
    /* Do not free these variables */
    msym_element_t *melements = NULL;
    msym_basis_function_t *mbfs = NULL;
    /* these are not mutable */
    const msym_symmetry_operation_t *msops = NULL;
    const msym_subgroup_t *msg = NULL;
    const msym_subrepresentation_space_t *msrs = NULL;
    const msym_character_table_t *mct = NULL;
    const msym_equivalence_set_t *mes = NULL;
    int mesl = 0;
    double *irrep = NULL;
    
    msym_basis_function_t *bfs = NULL;
    
    int msgl = 0, msopsl = 0, mlength = 0, msrsl = 0, mbfsl = 0, bfsl = 0;
        
    /* This function reads xyz files.
     * It initializes an array of msym_element_t to 0,
     * then sets the coordinates and name of the elements */
    int length = read_xyz(in_file, &elements);
    if(length <= 0) return -1;
    
    
    double (*psalcs)[bfsl] = NULL; // SALCs in matrix form, and input for symmetrization
    double *pcmem = NULL; // Some temporary memory
    int *pspecies = NULL;
    msym_partner_function_t *ppf = NULL;
    
    /* Create a context */
    msym_context ctx = msymCreateContext();
    
    if(NULL != thresholds){
        if(MSYM_SUCCESS != (ret = msymSetThresholds(ctx, thresholds))) goto err;
    }
    
    /* Use default thresholds otherwise call:
     * msymSetThresholds(msym_context ctx, msym_thresholds_t *thresholds); */
    
    /* Set elements */
    if(MSYM_SUCCESS != (ret = msymSetElements(ctx, length, elements))) goto err;
    
    /* Get elements msym elements */
    if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &mlength, &melements))) goto err;
    
    /* These are no longer needed, internal versions of these are kept in the context,
     * They are indexed in the same way that they have been allocated.
     * I.e. during orbital symmetrization or when getting the symmetrized LCAO,
     * the coefficients will correspond to the same indexing as "orbitals",
     * this is the main reason for the two levels of indirection */
    free(elements);  elements = NULL;
    
    /* Some trivial information */
    if(MSYM_SUCCESS != (ret = msymGetCenterOfMass(ctx,cm))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetRadius(ctx,&radius))) goto err;
    
//     printf("Molecule has center of mass [%lf; %lf; %lf] "
//            "and a radius of %lf\n",cm[0],cm[1],cm[2],radius);
//     
//     /* Find molecular symmetry */
    if(MSYM_SUCCESS != (ret = msymFindSymmetry(ctx))) goto err;
    
    /* Get the point group name */
    if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetSubgroups(ctx, &msgl, &msg))) goto err;
    
//     fprintf(stderr, "Found point group [0] %s with %d subgroups:\n",point_group, msgl);
//     fprintf(stderr, "\t [0] %s\n\t -------\n", point_group);
//     for(int i = 0; i < msgl;i++) fprintf(stderr,"\t [%d] %s\n",i+1,msg[i].name);
    int ssg = 0;
    
//     do {printf("\nChoose point group to use [0-%d]:",msgl);} while(scanf(" %d", &ssg) <= 0 || ssg < 0 || ssg > msgl);
//     if(ssg > 0){
//         ssg--;
//         printf("Selected point group %s\n",msg[ssg].name);
//         if(MSYM_SUCCESS != (ret = msymSelectSubgroup(ctx, &msg[ssg]))) goto err;
//         if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
//     }
    ssg = 0;
    
    char yn = 'n';
    
    /* Get elements msym elements */
    if(MSYM_SUCCESS != (ret = msymGetSymmetryOperations(ctx, &msopsl, &msops))) goto err;
    
    
    /* Set pointgroup to the C3v subgroup if it has XXX symmetry
     * using the same alignment as the original.
     * If specific axes are wanted the alignment axes/transform can be set using:
     * msym Get/Set Alignment Transform/Axes */
    if(0 == strncmp(point_group, "XXX", 3) && ssg == 0){
        //double transform[3][3];
        printf("Changing pointgroup from XXX -> C3v\n");
        //if(MSYM_SUCCESS != (ret = msymGetAlignmentTransform(ctx, transform))) goto err;
        if(MSYM_SUCCESS != (ret = msymSetPointGroupByType(ctx, MSYM_POINT_GROUP_TYPE_Cnv,3))) goto err;
        //if(MSYM_SUCCESS != (ret = msymSetAlignmentTransform(ctx, transform))) goto err;
        if(MSYM_SUCCESS != (ret = msymFindSymmetry(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
    }
    
    /* Retreive the symmetry operations */
    if(MSYM_SUCCESS != (ret = msymGetSymmetryOperations(ctx, &msopsl, &msops))) goto err;
        
    int symprint = 0;
        /* Symmetrize the molecule.
         * You can do this before orbital symmetrization as well,
         * but the permutations are already built, so you don't need to */
        if(MSYM_SUCCESS != (ret = msymSymmetrizeElements(ctx, &symerr))) goto err;
//         fprintf(stderr,"Molecule has been symmetrized to point group %s "
//                "with an error of %lf\n",point_group, symerr);
    
        if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &mlength, &melements))) goto err;
        if(mlength != length){ printf("Not possible!\n"); goto err;}
        
	printf("%d\n%s%15.8f\n",mlength, point_group,symerr);
        for(int i = 0;i < mlength;i++){
            printf("%s %12.9lf %12.9lf %12.9lf\n",
                   melements[i].name,
                   melements[i].v[0],
                   melements[i].v[1],
                   melements[i].v[2]);
        }
    
    msymReleaseContext(ctx);
    
    free(psalcs);
    free(pcmem);
    free(pspecies);
    free(ppf);
    free(bfs);
    free(elements);
    free(irrep);
    
    return ret;
err:
    free(psalcs);
    free(pcmem);
    free(pspecies);
    free(ppf);
    free(bfs);
    free(elements);
    free(irrep);
    error = msymErrorString(ret);
    fprintf(stderr,"Error %s: ",error);
    error = msymGetErrorDetails();
    fprintf(stderr,"%s\n",error);
    msymReleaseContext(ctx);
    return ret;
}

int main(int argc, const char * argv[]) {
    int ret = 1;
    if(argc == 2){
        ret = example(argv[1],NULL);
        fflush(stdout);
    } else {
        printf("usage msym_example <xyz-file>");
    }
    return ret;
}
