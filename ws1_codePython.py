import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

def fitPoisson(x,y):
	params=sm.GLM(y,x,family=sm.families.Poisson()).fit()
	pvalue=params.pvalues[0]
	fitted=params.predict(x)
	return fitted, pvalue

def tsvToArray(tsvPath):
	with open(tsvPath) as tsv:
		array=[line.strip().split() for line in tsv]
	return np.asarray(array)

def main(tsvPath):
	data=tsvToArray(tsvPath)
	data=data[1:len(data)]	
	samples=data[:,0]
	classes=data[:,1]
	hqcVals=np.asarray(data[:,2],dtype=np.int8)
	hqbVals=np.asarray(data[:,3],dtype=np.int8)
	hqtVals=np.asarray(data[:,4],dtype=np.int8)
	qtcVals=np.asarray(data[:,5],dtype=np.int8)
	qtbVals=np.asarray(data[:,6],dtype=np.int8)
	qttVals=np.asarray(data[:,7],dtype=np.int8)
	ages=np.asarray(data[:,8],dtype=np.float32)
	con=np.asarray(data[:,9],dtype=np.float32)
	black=(0,0,0)
	m_col=(72./255,118./255,255./255,120./255)
	m_col_opaque=(72./255,118./255,255./255)
	m_col_dark=(52./255,86./255,199./255,120./255)
	c_col=(255./255,99./255,71./255,120./255)
	c_col_opaque=(255./255,99./255,71./255)
	c_col_dark=(228./255,92./255,64./255,120./255)
	motherIndices=np.where(classes == "M")
	childIndices=np.where(classes == "C")	
	motherAges=ages[motherIndices]
	mother_hqt=hqtVals[motherIndices]
	childAges=con[childIndices]
	child_hqt=hqtVals[childIndices]
	motherFitted,motherP=fitPoisson(motherAges, mother_hqt)
	childFitted,childP=fitPoisson(childAges, child_hqt)

	ax1=plt.subplot2grid((10,10),(0,0),7,10,frameon=False)
	plt.plot(motherAges/365,motherFitted,linestyle="solid",lw=4,marker=None,color=m_col_opaque)
	plt.plot(childAges/365,childFitted,linestyle="solid",lw=4,marker=None,color=c_col_opaque)
	plt.plot(motherAges/365, mother_hqt,linestyle="None",marker="D",markerfacecolor=m_col,mec=m_col_dark, mew=2, markersize=10,label="mother p={p:.3f}".format(p=motherP))	
	plt.plot(childAges/365, child_hqt,linestyle="None",marker="o",markerfacecolor=c_col,mec=c_col_dark, mew=2, markersize=10, label="child p={p:.3f}".format(p=childP))
	plt.ylabel("number of point heteroplasmies",fontsize=12)
	plt.axis([15,60,-0.25,5.5],frameon=False)
	plt.axvline(x=15,ymin=0,ymax=5,color=black,lw=4)
	plt.axhline(y=-0.25,xmin=0,xmax=1,color=black,lw=6)	 
	plt.tick_params(direction="out", top=False, right=False, length=12,width=3, pad=10, labelsize=15)
	plt.legend(numpoints=1,loc=2,title="Poisson model",fontsize=12,handletextpad=0)
	
	ax2=plt.subplot2grid((10,10),(8,0),2,10,sharex=ax1,frameon=False)	
	plt.tick_params(top=False,bottom=False,right=False,left=False,labeltop=False,labelbottom=False,labelright=False,labelleft=False)
	plt.axis([15,60,0,1],frameon=False)
	plt.plot([min(childAges)/365.,max(childAges)/365.],[0.7,0.7],color=c_col_opaque,lw=2)
	plt.plot([min(motherAges)/365.,max(motherAges)/365.],[0.8,0.8],color=m_col_opaque,lw=2)
	plt.text(40,0.5,"collection",color=m_col_opaque,fontsize=14)
	plt.text(18,0.5,"fertilization",color=c_col_opaque,fontsize=14)
	plt.text(27,0.25,"maternal age (years)",color=black,fontsize=18)
	plt.show()

if __name__ == "__main__":
	main("ws1_data.txt")
