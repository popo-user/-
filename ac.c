#include<stdio.h>
#include<stdlib.h>
#include<time.h> 
#define LL long long
//#define int long long
int test1[10]={3,5,7,11};
int test2[10]={2,3,5,7,13,17,19};
struct big
{
	char s[1030];
	int len;
};
struct big32
{
	long long s[32];
}P;
typedef struct big  big;
typedef struct big32 big32;

void print(big x)
{
	printf("%d||",x.len);
	for(int i=x.len-1;~i;i--)
	{
		printf("%d",x.s[i]);
	}
	printf("\n");
}

void print32(big32 x)
{
	for(int i=31;~i;i--)
	{
		printf("%lld",x.s[i]);
	}
	printf("\n");
}
int up;
int ck(big x)
{
	big32 tx;
	for(int i=0;i<32;i++)
	{
		tx.s[i]=0;
		for(int j=0;j<32;j++)
		{
			if((i<<5|j)<x.len&&x.s[(i<<5)|j]) tx.s[i]|=(1ll<<j);
		}
	}
	for(int k=0;k<4;k++)
	{
		LL sum=0;
		for(int i=31;~i;i--)
		{
			sum=((sum<<32)+tx.s[i]%test1[k])%test1[k];
		}
		if(!sum) return 0; 
	}
	return 1;
}
static inline big mol (big x,big y)
{
	big32 tx,tz;
	for(int i=0;i<32;i++) tz.s[i]=0;
	for(int i=0;i<32;i++)
	{
		tx.s[i]=0;
		for(int j=0;j<32;j++)
		{
			if((i<<5|j)<x.len&&x.s[(i<<5)|j]) tx.s[i]|=(1ll<<j);
		}
	}
//	printf("%d\n",ty.s[0]&(1<<1));
	for(int i=y.len-1;~i;i--)
	{
		for(int j=31;~j;j--) 
		{
			tz.s[j]<<=1;
			if(j!=31&&(tz.s[j]&(1ll<<32)))
			{
				tz.s[j+1]|=1;
				tz.s[j]^=(1ll<<32);
			}
		}
//			printf("x");print32(tz);
		while(tz.s[31]&(1ll<<32))
		{
			for(int j=0;j<32;j++)
			{
				if(tz.s[j]<P.s[j]) {tz.s[j]+=(1ll<<32);tz.s[j+1]-=1;}
				tz.s[j]-=P.s[j];
			}
		}
			if(y.s[i])
			{
//				printf("%d ",i);
				for(int j=0;j<32;j++)
				{
					tz.s[j]+=tx.s[j];
					if(j!=31&&tz.s[j]&(1ll<<32))
					{
						tz.s[j]^=(1ll<<32);
						tz.s[j+1]++;
					}
				}
			}
			if(tz.s[31]&(1ll<<32))
			{
				for(int j=0;j<32;j++)
				{
					if(tz.s[j]<P.s[j]) {tz.s[j]+=(1ll<<32);tz.s[j+1]-=1;}
					tz.s[j]-=P.s[j];
				}
			}
	}
	int flow=0;
	for(int j=31;~j;j--)
	{
		if(tz.s[j]<P.s[j]) {flow=0;break;}
		else if(tz.s[j]>P.s[j]) {flow=1;break;}
	}
	if(flow)
	{
		for(int j=0;j<32;j++)
		{
			if(tz.s[j]<P.s[j]) {tz.s[j]+=(1ll<<32);tz.s[j+1]-=1;}
			tz.s[j]-=P.s[j];
		}
	}
	big z;
	for(int i=0;i<1026;i++) z.s[i]=0;
	for(int i=0;i<32;i++)
	{
		for(int j=0;j<32;j++)
		{
			if(tz.s[i]&(1<<j)) z.s[(i<<5)|j]=1;
		}
	}
	for(int i=1023;~i;i--)
	{
		if(z.s[i]) {z.len=i+1;break;}
	}
	
	return z;
}

big Pow_B(int low,big y)
{
	big x,res;int tot=0;
	
	for(int i=0;i<1025;i++) x.s[i]=0,res.s[i]=0;

	while(low)
	{
		x.s[tot++]=low&1;
		low>>=1;
	}
	x.len=tot;
	
	res.s[0]=1,res.len=1;
	
	up=0;	
	
	while(up<y.len)
	{
		if(y.s[up]&1) 
		{
			res=mol(res,x);		
		}
		x=mol(x,x);
		up++;
	}
	return res;

}

_Bool qus_B(big p)
{	
//	if(!p.s[0]) return 0;
	if(!ck(p)) return 0;
	big t=p;t.s[0]=0;big temp=t;
	
	
	int k = 0; while(!t.s[k]) k++;
	
	for(int i=0;i<t.len-k;i++) t.s[i]=t.s[i+k]; 	t.len-=k;
	
	for(int i=0;i<32;i++)
	{
		P.s[i]=0;
		for(int j=0;j<32;j++)
		{
			if(p.s[(i<<5)|j]) P.s[i]|=(1ll<<j);
		}
	}
	for(int i=0;i<5;i++)
	{
		big a=Pow_B(test2[i],t);
		big nxt=a;
		for(int j=0;j<k;j++)
		{
			nxt=mol(a,a);
			if((nxt.len==1&&nxt.s[0]==1))
			{
				if(!(a.len==1&&a.s[0]==1))//!=1
				{
					//!=-1 | temp
					if(a.len!=temp.len) return 0;
				
					for(int j=0;j<temp.len;j++)
					{
						if(a.s[i]!=temp.s[i]) return 0;
					}
				}	
			} 
			a=nxt;
		}
		if(!(a.len==1&&a.s[0]==1))return 0;
	}
	return 1;
}



big Rand_B()//1024
{
	big p;
	p.len=1024;
	p.s[1024]=0;p.s[1023]=1;//高位存高地址 
	for(int i=0;i<68;i++)//68*15=1020
	{
		int cur=rand();
		for(int j=0;j<15;j++)
		{
			p.s[i*15+j]=(cur&1)?1:0;
			cur>>=1;
		}
	}
	int cur=rand();
	for(int i=0;i<3;i++)
	{
		p.s[1020+i]=(cur&1)?1:0;
		cur>>=1;
	}
	p.s[0]=1;
	return p;
}


void work2()
{
	printf("a");
	int cnt=0;
	big p;
	while(cnt<1000)
	{
		cnt++;
		p=Rand_B(); 
		printf("%d ",cnt);
		if(qus_B(p)){print(p);break;}
	} 
	printf("%d",cnt);
}
int main()
{	
	
	srand((unsigned int)time(NULL));
	big p,t;
	for(int i=0;i<1025;i++) p.s[i]=0,t.s[i]=0; t.len=1024;p.len=1024;
	work2();
	
	return 0;	
}

