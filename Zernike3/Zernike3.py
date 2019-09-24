##     Zernike3.py was generated from Zernike.py using the 2to3 package. The original license is 
##     included below.

##     Zernike.py

##     This code takes a voxelised object and represents the shape of the object by computing a set
##     of Zernike moments.  It is also able to take a set of moments and from them reconstruct a
##     voxelised object.  Code is present that will load a pdb (protein databank) file and creat
##     a voxelised object from the protein structure.
##
##     The code was used in scientific research and as far as I know is correct.  However it is
##     mostly uncommented and can be improved in many ways.  For further details please see the
##     article:
##
##     Grandison S., Roberts C., Morris R. J. (2009)
##     The application of 3D zernike moments for the description of "model-free" molecular
##     structure, functional motion, and structural reliability
##     Journal of Computational Biology 16 (3) 487-500
##     DOI:10.1089/cmb.2008.0083
##
##     Or contact the author direct:  s.grandison@uea.ac.uk
##
##     Copyright (C) 2009  Scott Grandison, Richard Morris

##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.

##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.

##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.

import profile,random

import math,pickle

class StandardLibrary:
    "A class to contain some standard functions"
    def __init__(self):
        print("Created standard library")
    def CalcDistance(self,x1,y1,z1,x2,y2,z2):
        dist=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
        dist=dist**0.5
        return dist

class Moments:
    "A class to store Zernike moments"
    moments=[]
    invariants=[]
    def __init__(self):
        print("Created moments object")
    def ResizeMoments(self,size):
        NewMoments=[]
        for i in range(size+1):
            page=[]
            for j in range(size+1):
                line=[]
                for k in range(size+1):
                    line.append([0.0+0.0j,0.0+0.0j])
                page.append(line)
            NewMoments.append(page)
        if len(self.moments)>size+1:
            copyrange=size+1
        else:
            copyrange=len(self.moments)
        for i in range(copyrange):
            for j in range(copyrange):
                for k in range(copyrange):
                    NewMoments[i][j][k]=self.moments[i][j][k]
        self.moments=NewMoments
    def CalcInvariants(self):
        momrange=self.GetMomentSize()+1
        self.invariants=[]
        for (i,j) in [(a,b) for a in range(momrange) for b in range(a+1)]:
            momtot=0
            for k in range(-b,b+1,1):
                mom=self.GetMoment(i,j,k)
                if mom!="VOID":
                    momtot=mom.real**2+mom.imag**2
            self.invariants.append([i,j,momtot])
    def DisplayInvariants(self):
        for inv in self.invariants:
            if inv[2]!=0:
                print((inv[0],inv[1],inv[2]))
    def GetMomentSize(self):
        return len(self.moments)-1
    def GetMoment(self,n,l,m):
        val="VOID"
        index=0
        if m<0:
            index=1
        if ((n-l)%2 and 'Odd' or 'Even') =='Even' and l<=n and abs(m)<=l:
            val=self.moments[n][l][abs(m)][index]
        return val
    def SetMoment(self,n,l,m,val):
        if n>self.GetMomentSize():
            self.ResizeMoments(n)
        index=0
        if m<0:
            index=1
        self.moments[n][l][abs(m)][index]=val
    def DisplayAllMoments(self):
        momindex=len(self.moments)
        for i in range(momindex):
            for j in range(momindex):
                for k in range(-momindex,momindex):
                    val=self.GetMoment(i,j,k)
                    if val!="VOID" and val!=0:
                        print((i,j,k,val))
    def SaveMoments(self,filename):
        output=open(filename,"wb")
        pickle.dump(self.moments,output)
        output.close()
    def LoadMoments(self,filename):
        input=open(filename,"rb")
        self.moments=pickle.load(input)
        input.close()

class Voxels:
    "A class to store a voxelised grid"
    stan=StandardLibrary()
    voxels=[]
    infodictionary={}
    filename="NONE"
    filetype="NONE"
    filelines=[]
    newmaxextent=0.6 #Fraction of the grid to fit a structure into
    splatfrac=0.28 # Fraction of the grid to splat atoms over
    def __init__(self):
        print("Created a voxel object")
    def SetResolution(self,size):
        newvoxels=[]
        for i in range(size):
            page=[]
            for j in range(size):
                line=[]
                for k in range(size):
                    line.append(0.0)
                page.append(line)
            newvoxels.append(page)
        if len(self.voxels)>size:
            copyrange=size
        else:
            copyrange=len(self.voxels)
        for i in range(copyrange):
            for j in range(copyrange):
                for k in range(copyrange):
                    newvoxels[i][j][k]=self.voxels[i][j][k]
        self.voxels=newvoxels
    def GetResolution(self):
        return len(self.voxels)
    def Mirror(self):
        newvoxs=Voxels()
        res=self.GetResolution()
        newvoxs.SetResolution(res)
        for (i,j,k) in [(a,b,c) for a in range(res) for b in range(res) \
                            for c in range(res)]:
            newvoxs.SetVoxel(res-1-i,j,k,self.GetVoxel(i,j,k))
            #newvoxs.SetVoxel(res-1-i,j,k,0.0)
        return newvoxs
    def SetVoxel(self,x,y,z,value):
        if max(x+1,y+1,z+1)>self.GetResolution():
            self.SetResolution(max(x+1,y+1,z+1))
        self.voxels[x][y][z]=value
    def GetVoxel(self,x,y,z):
        return self.voxels[x][y][z]
    def DisplayAllVoxels(self):
        size=self.GetResolution()
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    val=self.GetVoxel(i,j,k)
                    if val!=0.0:
                        print((i,j,k,val))
    def LoadLines(self,filename):
        self.filename=filename
        f=open(filename,"r")
        lines=f.readlines()
        self.filename=filename
        self.filelines=lines
    def ConvertPDB(self):
        self.filetype="PDB"
        print("Loading PDB file")
        atoms=[]
        xtot=0;ytot=0;ztot=0;
        for line in self.filelines:
            if line[0:4]=="ATOM" or line[0:4]=="HETA":
                x=float(line[31:37])
                y=float(line[39:45])
                z=float(line[47:53])
                xtot+=x;ytot+=y;ztot+=z
                bfactor=float(line[61:66])/(8.0*math.pi**2)
                atom=line[77]
                electrons=8
                #Switch statement here
                atoms.append([x,y,z,bfactor,electrons])
            if line[0:5]=="TITLE":
                self.infodictionary["TITLE"]=line
            if line[0:6]=="AUTHOR":
                self.infodictionary["AUTHOR"]=line
        atomcount=len(atoms)
        print(("Detected ",atomcount," atoms"))
        xtot/=atomcount;ytot/=atomcount;ztot/=atomcount
        distancetable=[]
        for atom in atoms:
            atom[0]-=xtot
            atom[1]-=ytot
            atom[2]-=ztot
            dist=self.stan.CalcDistance(0,0,0,atom[0],atom[1],atom[2])
            distancetable.append(dist)
        distmax=max(distancetable)
        voxelsize=distmax/(self.GetResolution()*self.newmaxextent/2.0)
        vdw=1.7/voxelsize
        extmult=self.newmaxextent/distmax
        for atom in atoms:
            atom[0]*=extmult
            atom[1]*=extmult
            atom[2]*=extmult
        res=self.GetResolution()
        splatrange=int(self.splatfrac*self.GetResolution())
        for atom in atoms:
            xpos=0.5*(atom[0]+1)*res
            ypos=0.5*(atom[1]+1)*res
            zpos=0.5*(atom[2]+1)*res
            for xsplat in range(splatrange):
                for ysplat in range(splatrange):
                    for zsplat in range(splatrange):
                        xrel=-(splatrange-1)/2+xsplat
                        yrel=-(splatrange-1)/2+ysplat
                        zrel=-(splatrange-1)/2+zsplat
                        xvox=xrel+int(xpos)
                        yvox=yrel+int(ypos)
                        zvox=zrel+int(zpos)
                        dist=self.stan.CalcDistance(xvox,yvox,zvox,xpos,ypos,zpos)
                        if dist<float(splatrange*voxelsize/2):
                            sigma=atom[3]+vdw
                            val=atom[4]*math.exp(-dist**2/(2.0*(sigma)**2))
                            val/=(sigma*2*math.pi)**0.2
                            self.voxels[int(xvox)][int(yvox)][int(zvox)]+=val
        print(("Loaded "+self.filename))
    def Grid2DX(self,filename):
        print(("Writing to DX file:",filename))
        output=open(filename,"w")
        voxres=self.GetResolution()
        voxtot=voxres**3
        output.write("# Produced from ZClass\n")
        output.write("#\n#\n#\n")
        output.write("object 1 class gridpositions counts " \
                         +str(voxres)+" "+str(voxres)+" "+str(voxres)+"\n")
        output.write("origin 0.0e+00 0.0e+00 0.0e+0\n")
        output.write("delta 9.375000e-01 0.000000e+00 0.000000e+00\n")
        output.write("delta 0.000000e+00 9.375000e-01 0.000000e+00\n")
        output.write("delta 0.000000e+00 0.000000e+00 9.375000e-01\n")
        output.write("object 2 class gridconnections counts " \
                         +str(voxres)+" "+str(voxres)+" "+str(voxres)+"\n")
        output.write("object 3 class array type double rank 0 items " \
                         +str(voxtot)+" data follows\n")
        ct=0
        strg=""
        filledvoxels=0
        emptyvoxels=0
        for x in range(voxres):
            for y in range(voxres):
                for z in range(voxres):
                    strg+=str(self.voxels[x][y][z])+" "
                    if self.voxels[x][y][x]==0:
                        emptyvoxels+=1
                    else:
                        filledvoxels+=1
                    ct+=1
                    if ct==1:
                        ct=0
                        strg+="\n"
                    output.write(strg)
                    strg=""
        output.write('attribute "dep" string "positions"\n')
        output.write('object "regular positions regular connections" ' \
                         +'class field\n')
        output.write('component "positions" value 1\n')
        output.write('component "connections" value 2\n')
        output.write('component "data" value 3\n')
        print(("Filled voxels:",filledvoxels,"Empty voxels:",emptyvoxels))
        output.close()
    def DisplayInfo(self):
        print((self.infodictionary))
    def LoadStructure(self,filename):
        if filename[-3:]=="pdb":
            self.LoadLines(filename)
            self.ConvertPDB()
        else:
            print("I do not understand this filetype.")
    def SaveVoxels(self,filename):
        output=open(filename,"wb")
        pickle.dump([self.voxels,self.filelines,self.infodictionary],output)
        output.close()
    def LoadVoxels(self,filename):
        input=open(filename,"rb")
        tempload=pickle.load(input)
        self.voxels=tempload[0]
        self.filelines=tempload[1]
        self.infodictionary=tempload[2]
        input.close()

class Zernike:
    "A class for calculating Zernike moments from voxels"
    factdict={}
    bindict={}
    clmdict={}
    Qklvdict={}
    voxels=0
    axis=0
    axes=0
    tabSpaceSum=[]
    reconmoments=[]
    chidict={}
    def __init__(self):
        print("Created a Zernike calculator instance")
    def fac(self,n):
        try:
            val=self.factdict[n]
        except:
            a=lambda n:n-1 + abs(n-1) and a(n-1)*int(n) or 1
            val=a(n)
            self.factdict[n]=val
        return val
    def Binomial(self,n,k):
        try:
            val=self.bindict[(n,k)]
        except:
            val=self.fac(n)/(self.fac(n-k)*self.fac(k))
            self.bindict[(n,k)]=val
        return val
    def clm(self,l,m):
        try:
            val=self.clmdict[(l,m)]
        except:
            val=((2*l+1)*self.fac(l+m)*self.fac(l-m))**0.5/self.fac(l)
            self.clmdict[(l,m)]=val
        return val
    def Qklv(self,k,l,v):
        try:
            val=self.Qklvdict[(k,l,v)]
        except:
            coeff=(-1.0)**(k+v)/float(2**(2*k))*((2*l+4*k+3)/3.0)**0.5
            coeff*=self.Binomial(2*k,k)*self.Binomial(k,v)*self.Binomial(2*(k+l+v)+1,2*k)
            coeff/=self.Binomial(k+l+v,k)
            val=coeff
            self.Qklvdict[(k,l,v)]=val
        return val
    def sum6(self,n,l,m,nu,alpha,beta,u,mu,x,y,z):
        temp=0+0j
        for v in range(mu+1):
            r=2*(v+alpha)+u
            s=2*(mu-v+beta)+m-u
            t=2*(nu-alpha-beta-mu)+l-m
            if x=="BLAH":
                temp=temp+self.Binomial(mu,v)*self.tabSpaceSum[r][s][t]
            else:
                temp=temp+self.Binomial(mu,v)*self.axis[x][r+2]*self.axis[y][s+2]*self.axis[z][t+2]
        return temp
    def sum5(self,n,l,m,nu,alpha,beta,u,x,y,z):
        temp=0+0j
        for mu in range((l-m)/2+1):
            temp=temp+(-1)**mu*2**(-2*mu)*self.Binomial(l,mu)*self.Binomial(l-mu,m+mu)*self.sum6(n,l,m,nu,alpha,beta,u,mu,x,y,z)
        return temp
    def sum4(self,n,l,m,k,nu,alpha,beta,x,y,z):
        temp=0+0j
        for u in range(m+1):
            temp=temp+(-1)**(m-u)*self.Binomial(m,u)*(0+1j)**u*self.sum5(n,l,m,nu,alpha,beta,u,x,y,z)
        temp=temp.real-(temp.imag)*(0+1j)
        return temp
    def sum3(self,n,l,m,k,nu,alpha,x,y,z):
        temp=0+0j
        for beta in range(nu-alpha+1):
            temp=temp+self.Binomial(nu-alpha,beta)*self.sum4(n,l,m,k,nu,alpha,beta,x,y,z)
        return temp
    def sum2(self,n,l,m,k,nu,x,y,z):
        temp=0+0j
        for alpha in range(nu+1):
            temp=temp+self.Binomial(nu,alpha)*self.sum3(n,l,m,k,nu,alpha,x,y,z)
        return temp
    def sum1(self,n,l,m,k,x,y,z):
        k = int(round(k))
        temp=0+0j
        for nu in range(k+1):
            temp=temp+self.Qklv(k,l,nu)*self.sum2(n,l,m,k,nu,x,y,z)
        return temp
    def Chinlmrst(self,n,l,m,x=None,y=None,z=None):
        if x==None:
            k=(n-l)/2
            coeff=self.clm(l,m)*2.0**(-m)+0j
            coeff=coeff*0.75/3.1415926
            val=coeff*self.sum1(n,l,m,k,x,y,z)
        else:
            try:
                val=self.chidict[(n,l,m,x,y,z)]
            except:
                k=(n-l)/2
                coeff=self.clm(l,m)*2.0**(-m)+0j
                coeff=coeff*0.75/3.1415926
                val=coeff*self.sum1(n,l,m,k,x,y,z)
                self.chidict[(n,l,m,x,y,z)]=val
        return val
    def CalcRad(self,x,y,z):
        return (x**2+y**2+z**2)**0.5
    def CreateAxis(self,axis):
        res=self.voxels.GetResolution()
        self.axis=[]
        for i in range(res):
            temp=[i,-1+2*i/(float(res)),1.0]
            for j in range(self.order+1):
                temp.append((-1+2*i/(float(res)))**(j+1))
            self.axis.append(temp)
    def CreateAxes(self,axis,axes):
        self.axes=[]
        for xaxis in axis:
            for yaxis in axis:
                for zaxis in axis:
                    if self.CalcRad(xaxis[1],yaxis[1],zaxis[1])<=1.0:
                        self.axes.append([xaxis[0],yaxis[0],zaxis[0],xaxis[1],yaxis[1],zaxis[1]])
    def SpaceSum(self,a,b,g):
        tot=0.0
        diff=self.axis[4][1]-self.axis[3][1]
        for point in self.axes:
            x1=point[3];y1=point[4];z1=point[5];
            x2=x1+diff;y2=y1+diff;z2=z1+diff;
            temp=(x2**(a+1)-x1**(a+1))
            temp*=(y2**(b+1)-y1**(b+1))
            temp*=(z2**(g+1)-z1**(g+1))
            temp/=(a+1)*(b+1)*(g+1)
            temp*=self.voxels.GetVoxel(point[0],point[1],point[2])
            tot=tot+temp
        return tot
    def ClearSpace(self):
        pc=0
        indlist=[]
        totlen=len(self.axes)
        for ind in range(len(self.axes)):
            point=self.axes[ind]
            if self.voxels.GetVoxel(point[0],point[1],point[2])==0.0:
                pc=pc+1
                indlist.append(ind)
        indlist.reverse()
        if len(indlist)>0:
            del indlist[-1]
        for ind in indlist:
            del self.axes[ind]
        return [pc,totlen]
    def ConstructSpaceSumTable(self,order):
        for r in range(order+1):
            print(("Constructing order",r,"integral"))
            page=[]
            for s in range(order+1):
                line=[]
                for t in range(order+1):
                    val=0.0
                    if r+s+t<=order:
                        val=self.SpaceSum(r,s,t)
                    line.append(val)
                page.append(line)
            self.tabSpaceSum.append(page)
    def CalculateMoments(self,passedvoxels,order):
        self.voxels=passedvoxels
        self.order=order
        self.CreateAxis(self.axis)
        self.CreateAxes(self.axis,self.axes)
        clearstats=self.ClearSpace()
        self.ConstructSpaceSumTable(order)
        print(("Elements with zero:",clearstats[0],"  Total Elements:",clearstats[1]))
        print(("Removed:",clearstats[0]/float(clearstats[1]),"%"))
        calcedmoments=Moments()
        for n in range(order+1):
            print(("Doing:",n))
            for l in range(n+1):
                for m in range(l+1):
                    if ((n-l)%2 and 'Odd' or 'Even') =='Even':
                        chi=self.Chinlmrst(n,l,m)
                        negchi=0+0j
                        calcedmoments.SetMoment(n,l,m,chi)
                        if m>0:
                            negchi=(-1)**m*(chi.real-chi.imag*1j)
                            calcedmoments.SetMoment(n,l,-m,negchi)
        return calcedmoments
    def InitialiseReconstruction(self,passedmoments,order,resolution):
        print("Reconstructing voxel")
        self.order=order
        self.voxels=Voxels()
        self.voxels.SetResolution(resolution)
        self.reconmoments=passedmoments
        self.CreateAxis(self.axis)
        self.CreateAxes(self.axis,self.axes)
        self.reconvoxels=Voxels()
        self.reconvoxels.SetResolution(resolution)
        print("Done initialisation")
    def ReconstructVoxel(self,x,y,z):
        momorder=self.reconmoments.GetMomentSize()
        point=0.0
        for n in range(min(self.order,momorder)):
            for l in range(n+1):
                for m in range(l+1):
                    mom=self.reconmoments.GetMoment(n,l,m)
                    if mom!=0.0 and mom!="VOID":
                        point+=self.Chinlmrst(n,l,m,x,y,z)*mom
                    mom=self.reconmoments.GetMoment(n,l,-m)
                    if mom!=0.0 and mom!="VOID":
                        point+=self.Chinlmrst(n,l,-m,x,y,z)*mom
        return point.real
    def ReconstructAll(self):
        xold=10000
        for point in self.axes:
            xp=point[0];yp=point[1];zp=point[2]
            if xp!=xold:
                xold=xp
                print(("Doing: ",xold))
            self.reconvoxels.SetVoxel(xp,yp,zp,self.ReconstructVoxel(xp,yp,zp))
        return self.reconvoxels
    def InRegion(self,xyz):
        x=xyz[0]
        y=xyz[1]
        z=xyz[2]
        a=False
        currentvoxel=self.ReconstructVoxel(x,y,z)
        if currentvoxel>1.0:
            a=True
        return [a,currentvoxel]
    def FindAdjacent(self,xyz):
        x=xyz[0]
        y=xyz[1]
        z=xyz[2]
        adjlist=[]
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                for k in [-1,0,1]:
                    xd=i+x
                    yd=j+y
                    zd=k+z
                    coords=[xd,yd,zd]
                    try:
                        aaa=self.visited[str(coords)]
                    except:
                        if not(i==0 and j==0 and k==0):
                            adjlist.append([i+x,j+y,k+z])
        return adjlist
    def ResetChi(self):
        self.chidict={}
    def SaveChiDict(self,filename):
        output=open(filename,"wb")
        pickle.dump(self.chidict,output)
        output.close()
    def LoadChiDict(self,filename):
        input=open(filename,"rb")
        self.chidict=pickle.load(input)
        input.close()
    def ReconstructGrow(self):
        res=self.reconvoxels.GetResolution()
        halfres=int(float(res)/2.0)
        growthlist=[[halfres,halfres,halfres]]
        region=[]
        self.visited={}
        while (len(growthlist)>0):
            if len(region)%100==0:
                print(("Region size:",len(region),"Growth list:" \
                    ,len(growthlist),"Visited list:",len(self.visited)))
            toadd=growthlist[0]
            del(growthlist[0])
            region.append(toadd)
            prospective=self.FindAdjacent(toadd)
            for exam in prospective:
                self.visited[str(exam)]='VISITED'
                testelement=self.InRegion(exam)
                if testelement[0]==True:
                    growthlist.append(exam)
                    self.reconvoxels.SetVoxel(exam[0],exam[1],exam[2],testelement[1])
        print(("Region size:",len(region),"Growth list:",len(growthlist)))
        return self.reconvoxels

def TestAll():
    test=Moments()
    print((test.GetMomentSize()))
    print((test.moments))
    test.SetMoment(2,0,0,5+2j)
    test.DisplayAllMoments()
    test.SetMoment(3,1,-1,2+3j)
    test.DisplayAllMoments()
    test.SaveMoments("test.mom")
    newmoments=Moments()
    newmoments.LoadMoments("test.mom")
    print("Displaying loaded moments....")
    newmoments.DisplayAllMoments()
    print("\n\nTesting Voxel Stuff...")
    test2=Voxels()
    test2.SetResolution(1)
    print((test2.voxels))
    print((test2.GetResolution()))
    test2.SetVoxel(0,1,0,1.2)
    test2.SetVoxel(0,2,2,1.5)
    print((test2.voxels))
    test2.SetResolution(2)
    print("\n")
    print((test2.voxels))
    print(("Voxel val:",test2.GetVoxel(0,1,0)))
    test2.DisplayAllVoxels()
    test2.SetResolution(64)
    test2.LoadStructure("1FXZ.pdb")
    test2.Grid2DX("1FXZ.dx")
    test2.DisplayInfo()
    test2.SaveVoxels("1FXZ.vox")
    test2=[]
    newvox=Voxels()
    newvox.LoadVoxels("1FXZ.vox")
    newvox.DisplayInfo()
    test3=Zernike()
    willthiswork=test3.CalculateMoments(newvox,10)
    willthiswork.DisplayAllMoments()
    test3.InitialiseReconstruction(willthiswork,10,64)
    rvoxels=test3.ReconstructAll()
    rvoxels.Grid2DX("1FXZq.dx")

def TestSmall():
    newvox=Voxels()
    newvox.LoadVoxels("test.vox")
    newvox.DisplayInfo()
    test3=Zernike()
    willthiswork=test3.CalculateMoments(newvox,8)
    willthiswork.DisplayAllMoments()
    willthiswork.SaveMoments("testgrow.mom")
    test3.InitialiseReconstruction(willthiswork,4,64)
    rvoxels=test3.ReconstructGrow()
    rvoxels.Grid2DX("recon2.dx")

def TestGrow():
    willthiswork=Moments()
    willthiswork.LoadMoments("testgrow.mom")
    test3=Zernike()
    test3.InitialiseReconstruction(willthiswork,4,64)
    rvoxels=test3.ReconstructGrow()
    rvoxels=test3.ReconstructGrow()
    rvoxels.Grid2DX("recon2.dx")

#vox=Voxels()
#vox.SetResolution(64)
#vox.LoadStructure("1po5_protein_bindingSites_0HEM.pdb")
#vox.Grid2DX("test.dx")

#TestAll()
#TestSmall()
#TestGrow()
