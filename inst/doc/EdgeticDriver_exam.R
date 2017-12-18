### R code from vignette source 'EdgeticDriver_exam.Rnw'

###################################################
### code chunk number 1: EdgeticDriver
###################################################
library(EdgeticDriver)


###################################################
### code chunk number 2: EdgeticDriver
###################################################
#obtain the data for gene expression.
Input.exp=GetExampleData(exampleData="Exp.input")
# view first rows of data
head(Input.exp)
#obtain the sample label
Cancer_s=GetExampleData(exampleData="Cancer_s")
Normal_s=GetExampleData(exampleData="Normal_s")
#obtain the protein interaction network
Network=GetExampleData(exampleData="network")
#obtain the mutation data
Mut=GetExampleData(exampleData="mut")


###################################################
### code chunk number 3: EdgeticDriver
###################################################
#Identifying the dysregulated edges between cancer and normal samples
Input.exp=GetExampleData(exampleData="Exp.input")
# view first rows of data
head(Input.exp)
#obtain the sample label
Cancer_s=GetExampleData(exampleData="Cancer_s")
Normal_s=GetExampleData(exampleData="Normal_s")
#obtain the protein interaction network
Network=GetExampleData(exampleData="network")
#identify the dysregulated edges
DysCN=EdgeticDys_CN(Input.exp,Network,thr=0.01)


###################################################
### code chunk number 4: EdgeticDriver
###################################################
#Identifying the dysregulated edges among cancer samples
Input.exp=GetExampleData(exampleData="Exp.input")
# view first rows of data
head(Input.exp)
#obtain the sample label
Cancer_s=GetExampleData(exampleData="Cancer_s")
Normal_s=GetExampleData(exampleData="Normal_s")
#obtain the protein interaction network
Network=GetExampleData(exampleData="network")
#identify the dysregulated edges
Input.exp=Input.exp[,Cancer_s]
DysC=EdgeticDys_CN(Input.exp,Network,thr=0.01)


###################################################
### code chunk number 5: EdgeticDriver
###################################################
#Identifying the dysregulated edges in both two conditions
Input.exp=GetExampleData(exampleData="Exp.input")
# view first rows of data
head(Input.exp)
#obtain the sample label
Cancer_s=GetExampleData(exampleData="Cancer_s")
Normal_s=GetExampleData(exampleData="Normal_s")
#obtain the protein interaction network
Network=GetExampleData(exampleData="network")
#identify the dysregulated edges
DysCN=EdgeticDys_CN(Input.exp,Network,thr=0.01)
Input.exp2=Input.exp[,Cancer_s]
DysC=EdgeticDys_CN(Input.exp2,Network,thr=0.01)
Dys.net=EdgeticDys_both(DysCN,DysC)

###################################################
### code chunk number 6: EdgeticDriver
###################################################
#Identifying the mutation mediated-dysregulated edges in both two conditions
#obtain the data for gene expression.
Input.exp=GetExampleData(exampleData="Exp.input")
# view first rows of data
head(Input.exp)
#obtain the sample label
Cancer_s=GetExampleData(exampleData="Cancer_s")
Normal_s=GetExampleData(exampleData="Normal_s")
#obtain the protein interaction network
Network=GetExampleData(exampleData="network")
#obtain the mutation data
Mut=GetExampleData(exampleData="mut")
DysCN=EdgeticDys_CN(Input.exp,Network,thr=0.01)
Input.exp2=Input.exp[,Cancer_s]
DysC=EdgeticDys_CN(Input.exp2,Network,thr=0.01)
Dys.net=EdgeticDys_both(DysCN,DysC)
Driver.mut=getEdgeticDriver(Dys.net,Mut,alpha=0.05,n.sim=1000)

###################################################
### code chunk number 7: sessionInfo
###################################################
sessionInfo()
