#include<stdio.h>
#include<stdlib.h>
#include<time.h> 
#define LL long long
//#define int long long
int test1[30]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};//第一轮小素数筛 
int test2[10]={2,3,5,7,13,17,19};//检验用素数 
struct big
{
	char s[1030];
	int len;
};
struct big32
{
	long long s[32];
};
struct big60//压位 
{
	long long s[18];
}P;
typedef struct big  big;
typedef struct big32 big32;
typedef struct big60 big60;
int up;

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

void echo1(big x)
{
	for(int i=x.len-1;~i;i-=4)
	{
		int res=0;
		for(int j=0;j<4;j++)
		{
			res=(res<<1)+x.s[i-j];
		}
		printf("%X",res);
	}
	printf("\n");
}

void echo2(big x)//输出 
{
	FILE *fp;
	fp=fopen("ans.txt","w");
	for(int i=x.len-1;~i;i-=4)
	{
		int res=0;
		for(int j=0;j<4;j++)
		{
			res=(res<<1)+x.s[i-j];
		}
		printf("%X",res);
		fprintf(fp,"%X",res);
	}
	printf("\n"); 
	fclose(fp);
}
int ck(big x)//第一轮筛 
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
	for(int k=0;k<24;k++)
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
static inline big mol (big x,big y)//1024位乘法，取模 
{
	big60 tx,tz;
	for(int i=0;i<17;i++)
	{
		tx.s[i]=0;
		for(int j=0;j<60;j++)
		{
			if((i*60+j)<x.len&&x.s[i*60+j]) tx.s[i]|=(1ll<<j);
		}
	}
	tx.s[17]=0;
	for(int j=0;j<4;j++) 
	{
		if(1020+j<x.len&&x.s[1020+j]) tx.s[17]|=(1<<j);
	}
	for(int i=0;i<18;i++) tz.s[i]=0;
	
	for(int i=y.len-1;~i;i--)
	{
		for(int j=17;~j;j--) 
		{
			tz.s[j]<<=1;
			if(j!=17&&(tz.s[j]&(1ll<<60)))
			{
				tz.s[j+1]|=1;
				tz.s[j]^=(1ll<<60);
			}
		}
		while(tz.s[17]&(1<<4))
		{
			for(int j=0;j<18;j++)
			{
				if(tz.s[j]<P.s[j]) {tz.s[j]+=(1ll<<60);tz.s[j+1]-=1;}
				tz.s[j]-=P.s[j];
			}
		}
		if(y.s[i])
		{
			for(int j=0;j<18;j++)
			{
				tz.s[j]+=tx.s[j];
				if(j!=17&&tz.s[j]&(1ll<<60))
				{
					tz.s[j]^=(1ll<<60);
					tz.s[j+1]++;
				}
			}
		}
		if(tz.s[17]&(1<<4))
		{
			for(int j=0;j<18;j++)
			{
				if(tz.s[j]<P.s[j]) {tz.s[j]+=(1ll<<60);tz.s[j+1]-=1;}
				tz.s[j]-=P.s[j];
			}
		}
	}
	int flow=0;
	for(int j=17;~j;j--)
	{
		if(tz.s[j]<P.s[j]) {flow=0;break;}
		else if(tz.s[j]>P.s[j]) {flow=1;break;}
	}
	if(flow)
	{
		for(int j=0;j<18;j++)
		{
			if(tz.s[j]<P.s[j]) {tz.s[j]+=(1ll<<60);tz.s[j+1]-=1;}
			tz.s[j]-=P.s[j];
		}
	}
	big z;
	for(int i=0;i<1026;i++) z.s[i]=0;
	for(int i=0;i<17;i++)
	{
		for(int j=0;j<60;j++)
		{
			if(tz.s[i]&(1ll<<j)) z.s[i*60+j]=1;
		}
	}
	for(int j=0;j<4;j++)
	{
		if(tz.s[17]&(1<<j)) z.s[1020+j]=1;
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


_Bool qus_B(big p)//miller rabin
{	
	if(!ck(p)) return 0;
	
	big t=p;t.s[0]=0;big temp=t;//t=s-1 
	
	int k = 0; while(!t.s[k]) k++;
	
	for(int i=0;i<t.len-k;i++) t.s[i]=t.s[i+k]; 	t.len-=k;
	//s-1==t*(2^k)
	
	for(int i=0;i<17;i++)
	{
		P.s[i]=0;
		for(int j=0;j<60;j++)
		{
			if(p.s[i*60+j]) P.s[i]|=(1ll<<j);
		}
	}
	P.s[17]=0;
	for(int j=0;j<4;j++)
	{
		if(p.s[1020+j]) P.s[17]|=(1<<j);
	}
	//P为模数
	 
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


big Rand_B()//生成1024位大素数 
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
	p.s[0]=1;//必为奇数 
	return p;
}


void work2()
{
	clock_t start=clock();
	int cnt=0;
	big p;
	while(cnt<1000)
	{
		cnt++;
		p=Rand_B(); 
//		printf("%d ",cnt); 
		if(qus_B(p)){echo2(p);break;}
	} 
	clock_t end=clock();
	double dur=(double)(end-start)/CLOCKS_PER_SEC;
	printf("本次生成素数用时为%.2fs\n",dur);
}


int main()
{	
	srand((unsigned int)time(NULL));
	int op; 
	printf("请输入操作类型：(0.退出程序 1.测试给定大素数 2.产生一个随机大素数 均为1024位) "); 
	scanf("%d",&op);
	while(op)
	{
		if(op==1)
		{
			printf("程序将给出以下10个随机的1024位数，并判断其是否为素数\n"); 
			big p;
			for(int i=0;i<10;i++)
			{
				p=Rand_B();
				echo1(p);
				if(qus_B(p)) printf("此数为素数\n");
				else printf("此数不为素数\n");
			}
		}
		else
		{
			work2();
		}
		printf("请输入操作类型 ");
		scanf("%d",&op); 
	}
	
	return 0;	
}
