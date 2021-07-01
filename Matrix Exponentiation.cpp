struct Matrix
{
    vector<vll> a=vector<vll>(n,vll(n));
    Matrix operator *(Matrix other)
    {
        Matrix product;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                for(int k=0;k<n;k++)
                {
                    product.a[i][k]+=a[i][j]*other.a[j][k];
                }
            }
        }
        return product;
    }
};
 
Matrix power(Matrix a, ll k)
{
    Matrix product;
    for(int i=0;i<n;i++)
    {
        product.a[i][i]=1;
    }
    while(k>0)
    {
        if(k%2==1)
        {
            product=product*a;
        }
        a=a*a;
        k/=2;
    }
    return product;
}