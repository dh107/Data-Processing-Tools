/* --------------- Function File:Alloc memory for array ----------------- */
/* --------------- Include 1D and 2D -------------------------------------*/

void *alloc1d (size_t n1, size_t size)
{
        void *p;
        if ((p=malloc(n1*size)) == NULL)
                return NULL;
        return p;
}
float *alloc1dfloat (size_t n1)
{
        return (float*)alloc1d(n1, sizeof(float));
}
int *alloc1dint (size_t n1)
{
        return (int*)alloc1d(n1, sizeof(int));
}
void **alloc2d (size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;
    if((p=(void **)malloc(n2*n1*sizeof(void *)))==NULL)
        return NULL;
    if((p[0]=(void *)malloc(n2*n1*size))==NULL)
        return NULL;
    for(i2=0; i2<n2; i2++)
        p[i2]=(char *)p[0]+size*n1*i2;
    return p;
}
float **alloc2float (size_t n1, size_t n2)
{
    return (float **)alloc2d(n1, n2, sizeof(float));
}
void zero1float (float *p, int n1)
{
        int i;
        for (i=0; i<n1; i++)
                p[i] = 0;
}
void zero1int (int *p, int n1)
{
        int i;
        for (i=0; i<n1; i++)
                p[i] = 0;
}
void zero2float (float **p, int n1, int n2)
{
    int i,j;
    for(i=0; i<n2; i++)
    {
        for(j=0; j<n1; j++)
        {
            p[i][j] = 0;
        }
    }
}

