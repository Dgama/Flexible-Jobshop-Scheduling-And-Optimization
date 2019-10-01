from gurobipy import *
import numpy as np
import pandas as pd


#读取数据
products=pd.read_csv("Products_selected.csv")
multi=pd.read_csv("Multi_process.csv")

#索引
products_list=list(range(0,20))
MC_list=list(range(0,9))

#产品划分的array
array_products=np.asarray([products.loc[products.product_id==pid+1][['MC','duration']].values for pid in products_list])

#划分机器array
array_MC=np.asarray([products.loc[products.MC==mid+1][['product_id','process']].values for mid in MC_list])

#大M
bigM=1000000

#名称
name_list=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']

#创建模型
m=Model("PPC_project")

#完成时间变量定义
product_vars_list=[]
lower_bound=[]
for pid in products_list:
    for process_id in range(0,len(array_products[pid])):
        product_vars_list.append((pid,process_id))
        lower_bound.append(array_products[pid][process_id][1])
c=m.addVars(product_vars_list,lb=lower_bound,ub=24*365,vtype=GRB.INTEGER,name='c')

#工序限制
for pid in products_list:
    for process_id in range(0,len(array_products[pid])-1):
        m.addConstr(c[pid,process_id+1]-c[pid,process_id]>=array_products[pid][process_id+1][1])

#创建次序变量，包括和变量以及首位变量;次序变量个数由机器上加工的个数决定
sequence_vars_list=[]
for MC in MC_list:
    LOM=len(array_MC[MC])
    for i in range(0,LOM):
        for j in range(0,LOM):
            sequence_vars_list.append((MC,i,j))
sequence_vars=m.addVars(sequence_vars_list,vtype=GRB.BINARY,name='sequence_vars')

#单机器上约束
for MC in MC_list:
    for i in range(0,len(array_MC[MC])):
        m.addConstr(sequence_vars[MC,i,i]==0)
        # m.addConstr(sequence_vars.sum(MC, '*', i) - sequence_vars[MC, i, i]<=1)
        # m.addConstr(sequence_vars.sum(MC, i, '*') - sequence_vars[MC, i, i]<=1)

        #转化为时间约束
        ipid = array_MC[MC][i][0]-1
        iprocess=array_MC[MC][i][1]-1
        for j in range(0,len(array_MC[MC])):
            if i!=j:
                jpid = array_MC[MC][j][0]-1
                jprocess = array_MC[MC][j][1]-1
                m.addConstr(c[jpid,jprocess]-c[ipid,iprocess]+bigM*(1-sequence_vars[MC,i,j])>=array_products[jpid][jprocess][1])
    m.addConstr(sequence_vars.sum(MC,'*','*')==len(array_MC[MC]-1))

#特殊约束1-工序开始
m.addConstr(c[0,7]-array_products[0][7][1]==c[1,8])
m.addConstr(c[13,7]-array_products[13][7][1]==c[14,8])
m.addConstr(c[0,0]-array_products[0][0][1]==c[2,1]-array_products[2][1][1])
m.addConstr(c[13,0]-array_products[13][0][1]==c[15,0]-array_products[15][0][1])
m.addConstr(c[3,1]-array_products[3][1][1]==c[4,2]-array_products[4][2][1])
m.addConstr(c[19,2]-array_products[19][2][1]==c[18,7]-array_products[18][7][1])
m.addConstr(c[6,5]-array_products[6][5][1]==c[7,4]-48)
m.addConstr(c[11,0]-array_products[11][0][1]==c[12,5]-24)
m.addConstr(c[5,3]-array_products[5][3][1]==c[8,1]-array_products[8][1][1])
m.addConstr(c[9,2]-array_products[9][2][1]==c[10,7]-array_products[10][7][1])
m.addConstr(c[17,6]-array_products[17][6][1]==c[18,8])
m.addConstr(c[0,len(array_products[0])-1]<=c[13,len(array_products[13])-1]-1)

#特殊约束2-开始、完工时间约束
m.addConstr(c[13,0]-array_products[13][0][1]>=24*14)
m.addConstr(c[14,0]-array_products[14][0][1]>=24*14)
m.addConstr(c[0,len(array_products[0])-1]-array_products[0][len(array_products[0])-1][1]<=24*60)
m.addConstr(c[3,len(array_products[3])-1]-array_products[3][len(array_products[3])-1][1]<=24*91)
m.addConstr(c[4,len(array_products[4])-1]-array_products[4][len(array_products[4])-1][1]<=24*152)
m.addConstr(c[19,0]-array_products[19][0][1]>=24*120)
m.addConstr(c[19,len(array_products[19])-1]-array_products[19][len(array_products[19])-1][1]<=24*181)
m.addConstr(c[18,len(array_products[18])-1]-array_products[18][len(array_products[18])-1][1]<=24*334)
m.addConstr(c[6,0]-array_products[6][0][1]>=24*151)
m.addConstr(c[6,len(array_products[6])-1]-array_products[6][len(array_products[6])-1][1]<=24*201)
m.addConstr(c[7,0]-array_products[7][0][1]>=24*59)
m.addConstr(c[7,len(array_products[7])-1]-array_products[7][len(array_products[7])-1][1]<=24*242)
m.addConstr(c[9,0]-array_products[9][0][1]>=24*151)
m.addConstr(c[9,len(array_products[9])-1]-array_products[9][len(array_products[9])-1][1]<=24*334)
m.addConstr(c[10,len(array_products[10])-1]-array_products[10][len(array_products[10])-1][1]<=24*334)
m.addConstr(c[11,len(array_products[11])-1]-array_products[11][len(array_products[11])-1][1]<=24*242)
m.addConstr(c[12,len(array_products[12])-1]-array_products[12][len(array_products[12])-1][1]<=24*242)
m.addConstr(c[5,0]-array_products[5][0][1]>=24*78)
m.addConstr(c[5,len(array_products[5])-1]-array_products[5][len(array_products[5])-1][1]<=24*303)
m.addConstr(c[8,0]-array_products[8][0][1]>=24*78)
m.addConstr(c[8,len(array_products[8])-1]-array_products[8][len(array_products[8])-1][1]<=24*303)
m.addConstr(c[17,0]-array_products[17][0][1]>=24*31)
m.addConstr(c[17,len(array_products[17])-1]-array_products[17][len(array_products[17])-1][1]<=24*273)
m.addConstr(c[2,len(array_products[2])-1]-array_products[2][len(array_products[2])-1][1]<=24*335)
m.addConstr(c[15,len(array_products[15])-1]-array_products[15][len(array_products[15])-1][1]<=24*335)
m.addConstr(c[16,len(array_products[16])-1]-array_products[16][len(array_products[16])-1][1]<=24*335)

#求最长完成时间
product_gurobi_vars_list=[]
for i in range(0,len(product_vars_list)):
    product_gurobi_vars_list.append(c[product_vars_list[i][0],product_vars_list[i][1]])
max_time=m.addVar(lb=0,name='max_time')
m.addConstr(max_time==max_(product_gurobi_vars_list))

#求解
m.setObjective(max_time,GRB.MINIMIZE)
m.optimize()

#输出排序
for i in product_vars_list:
    complete_time=m.getVarByName('c['+str(i[0])+','+str(i[1])+']').x
    begin_time=complete_time-array_products[i[0]][i[1]][1]
    begin_day=int(begin_time/24)
    begin_hour=begin_time-24*begin_day
    print('Product',name_list[i[0]],'Process',i[1],'Begins at day',begin_day,'hour',begin_hour)

total_days=int(m.ObjVal/24)
left_hours=m.ObjVal-total_days*24
print('Total Time is :',total_days,'days and',left_hours,'hours')