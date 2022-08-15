#include<stdio.h>
#include<math.h>
#include<stdlib.h>
void gauss_seidel();
void gauss_elim();
void main()
{
    int c;
    printf("\t\tEnter the Gauss method for solving linear equations\n");
    printf("\t\t\t 1.Gauss-Seidel Iteration Method\n");
    printf("\t\t\t 2.Gauss elimination method\n");
    printf("\t\t");
    scanf("%d",&c);
    switch(c)
    {
        case 1:gauss_seidel();
               break;
        case 2:gauss_elim();
               break;
        default:printf("\t\tInvalid choice");
    }

}
void gauss_seidel()
{
	float a[10][10],b[10],x[10],xn[10],epp=0.00001,sum;
	int i,j,n,flag,f,k;
	printf("Enter number of variables: ");
	scanf("%d",&n);
	printf("Enter the coefficients row-wise: ");
	for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("x[%d][%d]: ",(i+1),(j+1));
            scanf("%f",&a[i][j]);
        }
    }
    //checking for row dominance
    flag=0;
	for(i=0;i<n;i++)
	{
		sum=0;
		for(j=0;j<n;j++)
        {
            if(i!=j)
				sum+=fabs(a[i][j]);
			if(sum>fabs(a[i][i]))
				flag=1;
        }
	}
	//checking for column dominance
	if(flag==1)
		flag=0;
	for(j=1;j<n;j++)
	{
		sum=0;
		for(i=1;i<n;i++)
        {
            if(i!=j)
				sum+=fabs(a[i][j]);
			if(sum>fabs(a[j][j]))
				flag=1;
        }
	}

	if(flag==1)
	{
		printf("The coefficient matrix is not diagonally dominant \n");
		printf("The Gauss-Jacobi method does not converge surely");
		exit(0);
	}
	printf("Enter right hand constant values: ");
	for(i=0;i<n;i++)
    {
        printf("b[%d]: ",(i+1));
        scanf("%f",&b[i]);
    }
	for(i=0;i<n;i++)
		x[i]=0; //initialize
    sum=0;
	do{
		for(i=0;i<n;i++)
		{
			sum=b[i];
			for(j=0;j<n;j++)
			{
				if(j<i)
					sum-=a[i][j]*xn[j];
				else if(j>i)
					sum-=a[i][j]*x[j];
                xn[i]=sum/a[i][i];
			}
		}
		flag=0; // indicates |x[i]-xn[i]|<epp for all i
		for(i=0;i<n;i++)
        {
			if(fabs(x[i]-xn[i])>epp)
				flag=1;
        }
		if(flag==1)
        {
            for(i=0;i<n;i++)
				x[i]=xn[i]; // reset x[i]
        }
	}while(flag==1);
	printf("Solution is \n");
	for(i=0;i<n;i++)
		printf("x%d=%8.5f\t",(i+1),xn[i]);
}
void gauss_elim()
{
    int i,j,k,n;
    float A[20][20],c,x[10],sum=0.0;
    printf("\nEnter the order of matrix: ");
    scanf("%d",&n);
    printf("\nEnter the elements of augmented matrix row-wise:\n\n");
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=(n+1); j++)
        {
            printf("A[%d][%d] : ", i,j);
            scanf("%f",&A[i][j]);
        }
    }
    for(j=1; j<=n; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i=1; i<=n; i++)
        {
            if(i>j)
            {
                c=A[i][j]/A[j][j];
                for(k=1; k<=n+1; k++)
                {
                    A[i][k]=A[i][k]-c*A[j][k];
                }
            }
        }
    }
    x[n]=A[n][n+1]/A[n][n];
    /* this loop is for backward substitution*/
    for(i=n-1; i>=1; i--)
    {
        sum=0;
        for(j=i+1; j<=n; j++)
        {
            sum=sum+A[i][j]*x[j];
        }
        x[i]=(A[i][n+1]-sum)/A[i][i];
    }
    printf("\nThe solution is: \n");
    for(i=1; i<=n; i++)
    {
        printf("\nx%d=%f\t",i,x[i]); /* x1, x2, x3 are the required solutions*/
    }
}


