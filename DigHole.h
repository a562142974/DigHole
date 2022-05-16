#include<fstream>
#include<iostream>
#include<iomanip>
#include<string>
#include<vector>
#include"math.h"
#include <algorithm>

//点
class point{
    public:
    double x,y;
    int index;


    int *Eindexes;//用于存储含有这个点的单元的索引
    int nEs=0;//含有这个点的单元个数


    point* next=nullptr;
    //默认构造函数
    point(){Eindexes=new int[20];nEs=0;};
    //构造函数
    point(double nx,double ny,int nindex);
    //赋值函数
    void set(double nx,double ny,int nindex);
    //直接读
    void set(std::ifstream &fin,int index);
    //存储单元索引
    void fillEindexes(int nEindex);
};

point::point(double nx,double ny,int nindex):x(nx),y(ny),index(nindex){
    Eindexes=new int[10];
    nEs=0;
}

void point::set(double nx,double ny,int nindex){
    x=nx;
    y=ny;
    index=nindex;
}

void point::set(std::ifstream &fin,int nindex){
    fin>>x>>y;
    index=nindex;
};

void point::fillEindexes(int nEindex){
    Eindexes[nEs]=nEindex;
    nEs+=1;
}

//背景网格
class box{
    public:
    point* next=nullptr;
    box()=default;
};


class DigHole{
    public:

    //所有主网格点
    point** APts;
    //所有单元
    int** Elmts;

    //点数
    int N_pts;
    //单元数
    int N_elmts;

    //object
    //object点数
    int N_pts_obj;
    double** APts_obj;
    double r,xc,yc;//obj的外接圆半径，圆心坐标    
    std::ifstream object;
    

    //所有盒子
    box** boxes;
    int Nbx,Nby;

    //关于整个domain
    double xl,xr,yd,yu;


    //构造函数
    DigHole(std::string filename_maingrids,std::string filename_object);

    //挖洞
    void Dig(std::string filename_obj,std::string filename_output);//由文件输入新的obj位置挖洞
    void Dig(double** nAPts_obj,std::string filename_output);//由数组输入新的obj位置挖洞
    void Dig(std::string filename_output);//根据现有变量挖洞

    private:
    //不同的具体挖洞算法
    void Raydiscriminant(std::vector<point*> &neis,std::vector<int> &index_pts_dlted,std::vector<int> &index_elmts_dlted);
};

DigHole::DigHole(std::string filename_maingrids,std::string filename_object){
    std::ifstream maingrids(filename_maingrids);
    object.open(filename_object);
    //读取主网格点数和单元数
    std::string firstline;
    getline(maingrids,firstline);
    int starpos=0,length=0,k=0;
    while (firstline[k]<48 || firstline[k]>57)
    {
        k+=1;
    }
    starpos=k;
    while (firstline[k]>=48 && firstline[k]<=57)
    {   
        length+=1;
        k+=1;
    }
    N_pts=std::stoi(firstline.substr(starpos,length));
    length=0;
    while (firstline[k]<48 || firstline[k]>57)
    {
        k+=1;
    }
    starpos=k;
        while (firstline[k]>=48 && firstline[k]<=57)
    {   
        length+=1;
        k+=1;
    }
    N_elmts=std::stoi(firstline.substr(starpos,length));

    //APts赋值
    APts=new point*[N_pts];
    for (int i = 0; i < N_pts; i++)
    {
        APts[i]=new point;
        APts[i]->set(maingrids,i+1);
    }
    
    //Elmts赋值
    Elmts=new int*[N_elmts];
    for (int i = 0; i < N_elmts; i++)
    {
        Elmts[i]=new int[3];
        maingrids>>Elmts[i][0];
        maingrids>>Elmts[i][1];
        maingrids>>Elmts[i][2];

        //赋值点的索引
        for (int j = 0; j < 3; j++)
        {
            APts[Elmts[i][j]-1]->fillEindexes(i+1);               
        }
    }

    //读取obj外延点数
    length=0;k=0;starpos=0;
    getline(object,firstline);
    while (firstline[k]<48 || firstline[k]>57)
    {
        k+=1;
    }
    starpos=k;
    while (firstline[k]>=48 && firstline[k]<=57)
    {   
        length+=1;
        k+=1;
    }
    N_pts_obj=std::stoi(firstline.substr(starpos,length));

    //APts_object;
    APts_obj=new double*[N_pts_obj];
    for (int i = 0; i < N_pts_obj; i++)
    {
        APts_obj[i]=new double[2];
        object>>APts_obj[i][0];
        object>>APts_obj[i][1];
    }
    
    //将所有点装盒
    //先计算xl,xr,yd,yu,假设APts的前四个点是边界点

    xl=APts[N_pts-4]->x;xr=APts[N_pts-3]->x;yu=APts[N_pts-1]->y;yd=APts[N_pts-4]->y;


    //计算r
    r=-1;
    for (int i = 0; i <=N_pts_obj-2; i++)
    {
        for (int k = i+1; k < N_pts_obj; k++)
        {
            r=sqrt(pow(APts_obj[i][0]-APts_obj[k][0],2)+pow(APts_obj[i][1]-APts_obj[k][1],2))>r ? sqrt(pow(APts_obj[i][0]-APts_obj[k][0],2)+pow(APts_obj[i][1]-APts_obj[k][1],2)):r;
        }
    }

    //装盒
    //初始化盒
    Nbx=(xr-xl)/r+1;
    Nby=(yu-yd)/r+1;
    boxes=new box*[Nbx];
    for (int i = 0; i < Nbx; i++)
    {
        boxes[i]=new box[Nby];
    }

    //装
    int m,n;
    for (int i = 0; i < N_pts; i++)
    {
        m=(APts[i]->x-xl)/r;n=(APts[i]->y-yd)/r;
        APts[i]->next=boxes[m][n].next;
        boxes[m][n].next=APts[i];
    }
}


void DigHole::Dig(std::string filename_output){
    //确定obj在哪些盒子
    int il=(APts_obj[0][0]-xl)/r,ir=(APts_obj[0][0]-xl)/r,jd=(APts_obj[0][1]-yd)/r,ju=(APts_obj[0][1]-yd)/r;
    for (int i = 1; i < N_pts_obj; i++)
    {
        il=(APts_obj[i][0]-xl)/r<il ? (APts_obj[i][0]-xl)/r:il;
        ir=(APts_obj[i][0]-xl)/r>ir ? (APts_obj[i][0]-xl)/r:ir;
        jd=(APts_obj[i][1]-xl)/r<jd ? (APts_obj[i][1]-xl)/r:jd;
        ju=(APts_obj[i][1]-xl)/r>ju ? (APts_obj[i][1]-xl)/r:ju;
    }

    //将盒子中的点压入neis中
    point* p;
    std::vector<point*> neis;
    for (int i = il; i <= ir; i++)
    {
        for (int j = jd; j <= ju; j++)
        {
            p=boxes[i][j].next;
            while(p!=nullptr){
                neis.push_back(p);
                p=p->next;
            }
        }
    }

    //调用具体dig算法
    std::vector<int> index_pts_dlted;
    std::vector<int> index_elmts_dlted;
    std::cout<<"开始挖洞..."<<std::endl;
    Raydiscriminant(neis,index_pts_dlted,index_elmts_dlted);
    std::cout<<"挖洞结束..."<<std::endl;

    //根据std::vector<int> &index_pts_dlted,std::vector<int> &index_elmts_dlted输出plt文件
    int* flag_pts_dlted;int* flag_elmts_dlted;
    flag_pts_dlted=new int[N_pts]();
    flag_elmts_dlted=new int[N_elmts]();
    for (int i = 0; i < index_pts_dlted.size(); i++)
    {
        flag_pts_dlted[index_pts_dlted[i]-1]=1;
    }
    for (int i = 0; i < index_elmts_dlted.size(); i++)
    {
        flag_elmts_dlted[index_elmts_dlted[i]-1]=1;
    }
    int N_pts_remained=N_pts-index_pts_dlted.size();
    int N_elmts_remained=N_elmts-index_elmts_dlted.size();


    //开始输出
    std::cout<<"正在输出到文件..."<<std::endl;
    std::ofstream fout(filename_output);
    //fout<<"ZONE N="+std::to_string(N_pts_remained)+", E="+std::to_string(N_elmts_remained)+", F=FEPOINT, ET=TRIANGLE"<<std::endl;
    fout<<"ZONE N="+std::to_string(N_pts)+", E="+std::to_string(N_elmts)+", F=FEPOINT, ET=TRIANGLE"<<std::endl;
    for (int i = 0; i < N_pts; i++)
    {
        //if(flag_pts_dlted[i]==1) continue;
        fout<<APts[i]->x<<" "<<APts[i]->y<<std::endl;
    }
    for (int i = 0; i < N_elmts; i++)
    {
        if(flag_elmts_dlted[i]==1) {
            fout<<1<<" "<<1<<" "<<1<<std::endl;
            //continue;
        }
        else{
        fout<<Elmts[i][0]<<" "<<Elmts[i][1]<<" "<<Elmts[i][2]<<std::endl;}
    }


    std::string temp;
    object.seekg(0);
    while (getline(object,temp))
    {
        fout<<temp<<std::endl;
    }
    fout.close();
    std::cout<<"完成"<<std::endl;
};



void DigHole::Raydiscriminant(std::vector<point*> &neis,std::vector<int> &index_pts_dlted,std::vector<int> &index_elmts_dlted){
    //in:std::vector<point*> &neis
    //out:std::vector<int> &index_pts_dlted,std::vector<int> &index_elmts_dlted
    //默认objct中的点按照顺时针或者逆时针排列
    //找index_pts_dlted
    int n_intersec=0;
    double x1,y1,x2,y2,x_intersec,y_intersec,x_min,x_max,y_min,y_max;
    double k=1;
    for (int i = 0; i < neis.size(); i++)
    {
        for (int j = 0; j < N_pts_obj; j++)
        {
            x1=APts_obj[j][0];y1=APts_obj[j][1];
            if(j!=N_pts_obj-1){
                x2=APts_obj[j+1][0];y2=APts_obj[j+1][1];
            }else{
                x2=APts_obj[0][0];y2=APts_obj[0][1];
            }

            x_min= x1<x2 ? x1:x2;
            x_max=x1>x2 ? x1:x2;
            y_min=y1<y2 ? y1:y2;
            y_max=y1>y2 ? y1:y2;
            x_intersec=(k*neis[i]->x-neis[i]->y-(y2-y1)/(x2-x1)*x1+y1)/(k-(y2-y1)/(x2-x1));
            y_intersec=1/((x2-x1)/(y2-y1)-1/k)*(neis[i]->x-neis[i]->y+y1/(y2-y1)*(x2-x1)-x1);
            if((x_intersec<x_max && x_intersec>x_min && x_intersec>neis[i]->x) | \
                (y_intersec<y_max && y_intersec>y_min && y_intersec>neis[i]->y)) \
            n_intersec+=1;
        }
        //判断n_intersec奇偶性
        if(n_intersec/2*2!=n_intersec) index_pts_dlted.push_back(neis[i]->index);
        n_intersec=0;
    }


    //根据index_pts_dlted找index_elmts_dlted

    
    //此处找点所在单元索引会导致单元索引重复，但不会影响后续输出，此处可以更改，根据vector语法
    int index;
    for (int i = 0; i < index_pts_dlted.size(); i++)
    {
        for (int j = 0; j < APts[index_pts_dlted[i]-1]->nEs; j++)
        {
            index=APts[index_pts_dlted[i]-1]->Eindexes[j];
            if(std::find(index_elmts_dlted.begin(),index_elmts_dlted.end(),index)==index_elmts_dlted.end()){
            index_elmts_dlted.push_back(index);}
        }
    }
}