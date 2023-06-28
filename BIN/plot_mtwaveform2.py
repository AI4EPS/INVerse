# plot_mtwaveform.py
# Script to plot waveform fit
# Yuexin Li and Gabriel Rogow-Patt
# First version: April, 2019
# Second version: 2021

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
#from obspy.imaging.beachball import beach,beachball
import math
import fpsol
import pandas as pd
#import sys
#sys.path.append('./')
from mopad_wrapper import beach
from numpy import linalg as la
import argparse
from matplotlib.backends.backend_pdf import PdfPages


parser = argparse.ArgumentParser()
parser.add_argument('--nodisplaywindow', action='store_true', help='prevents the plot(s) from opening any display window(s)')
parser.add_argument('--nosavefig', action='store_true', help='prevents the plot(s) from saving to any file(s)')
parser.add_argument('--noplotazimuth', action='store_true', help='prevents the station names from plotting on the beachball(s)')
parser.add_argument('--plotseparatebeachball', action='store_true', help='saves extra file(s) with the beachball diagram(s)')
parser.add_argument('--plottape', action='store_true', help='saves extra file(s) with the Tape Plot diagram(s)')
parser.add_argument("mtinv_file", action='extend', nargs='+')
args = parser.parse_args()

# Python script to substitute "tdmt_plot_gmt5.perl mtinvout"


moscale=1    #Units: Nm                                     #scale value for moment estimates and Mw

# Plotting while reading
for mtinput in args.mtinv_file:
    fig=plt.figure(figsize=(9,6))
    mark=0
    index=-1
    ncol=4
    DT=[];DR=[];DZ=[] #List for data
    ST=[];SR=[];SZ=[]
    filenameazlist=[]
    f=open(mtinput,'r')
    linelist = f.readlines()
    if len(linelist[-1])!=0:
        linelist.append('')
    #for line in f.readlines():
    with PdfPages('fm_waveform(' + mtinput + ').pdf') as pdf:
        for line in linelist:
            line=line.strip()
            if len(line)==0:
                mark=0
                index=index+1
                if index>0:
                    index=index-1
                    nstascale=nsta/4
                    ax3 = fig.add_subplot(nsta,ncol,1+(index%8)*ncol)
                    trans = transforms.blended_transform_factory(ax3.transData, ax3.transAxes)
                    #plt.subplot(nsta,4,1+(index-1)*4)
                    plt.plot(DT,'k',linewidth=.75)
                    #plt.plot(ST,'r--',linewidth=.75, scaley=False)
                    plt.axis('off')
                    plt.text(-.1,.5,filename,fontsize=8,transform=ax3.transAxes,horizontalalignment='right')
                    if (index%8)==0:
                        plt.title('Tangential')
                    text_str='Distance = %d km  Azimuth = %d  Max Amp = %.2e cm  Zcorr = %d  VR = %d' % (dist, Az, np.max(np.abs([DT,DR,DZ])), Zcor, VR)
                    plt.text(0,-.15, text_str,fontsize=8,transform=trans)


                    ax6 = fig.add_subplot(nsta,ncol,2+(index%8)*ncol,sharey=ax3)
                    #plt.subplot(nsta,ncol,2+(index-1)*ncol)
                    plt.plot(DR,'k',linewidth=.75)
                    #plt.plot(SR,'r--',linewidth=.75, scaley=False)
                    plt.axis('off')
                    if (index%8)==0:
                        plt.title('Radial')


                    ax2 = fig.add_subplot(nsta,ncol,3+(index%8)*ncol,sharey=ax3)
                    #plt.subplot(nsta,ncol,3+(index-1)*ncol)
                    plt.plot(DZ,'k',linewidth=.75)
                    #plt.plot(SZ,'r--',linewidth=.75, scaley=False)
                    plt.axis('off')
                    if (index%8)==0:
                        plt.title('Vertical')

                    trans2 = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
                    if DZ.index(min(DZ))<len(DZ)/2:
                        plt.text(npts-.125*npts,.3-.1*nstascale,str(round(.25*npts*dt,2)) + ' sec',fontsize=5,transform=trans2,horizontalalignment='center')
                        plt.plot([npts-.25*npts,npts],[.3,.3],'k',transform=trans2,linewidth=.75)
                        plt.plot([npts-.25*npts,npts-.25*npts],[.28,.32],'k',transform=trans2,linewidth=.75)
                        plt.plot([npts,npts],[.28,.32],'k',transform=trans2,linewidth=.75)
                        plt.axis('off')
                    else:
                        plt.text(.125*npts,.3-.1*nstascale,str(round(.25*npts*dt,2)) + ' sec',fontsize=5,transform=trans2,horizontalalignment='center')
                        plt.plot([.25*npts,0],[.3,.3],'k',transform=trans2,linewidth=.75)
                        plt.plot([.25*npts,.25*npts],[.28,.32],'k',transform=trans2,linewidth=.75)
                        plt.plot([0,0],[.28,.32],'k',transform=trans2,linewidth=.75)
                        plt.axis('off')

                    if (index%8)==0:
                        if nsta==2 and nsta1==False or nsta==3:
                            ax4 = fig.add_subplot(nsta,ncol,ncol)
                            plt.text(.2,1.1,'Depth = ' + str(depth),transform=ax4.transAxes)
                            plt.text(.2,1.1-.2*nstascale,'Strike = ' + str(round(st1)) + ' ; ' + str(round(st2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-.4*nstascale,'Rake = ' + str(round(rk1)) + ' ; ' + str(round(rk2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-.6*nstascale,'Dip = ' + str(round(dp1)) + ' ; ' + str(round(dp2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-.8*nstascale,'Mo = {:0.2e}'.format(Motot),transform=ax4.transAxes)
                            plt.text(.2,1.1-1.0*nstascale,'Mw = ' + str(round(Mw,2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-1.2*nstascale,'Percent DC = ' + str(round(pdc,2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-1.4*nstascale,'Percent CLVD = ' + str(round(pclvd,2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-1.6*nstascale,'Percent ISO = ' + str(round(piso,2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-1.8*nstascale,'Var. Red. = ' + str(round(VRtot,2)),transform=ax4.transAxes)
                            plt.text(.2,1.1-2.0*nstascale,'Mxx = {:0.3e}'.format(mxx),transform=ax4.transAxes)
                            plt.text(.2,1.1-2.2*nstascale,'Mxy = {:0.3e}'.format(mxy),transform=ax4.transAxes)
                            plt.text(.2,1.1-2.4*nstascale,'Mxz = {:0.3e}'.format(mxz),transform=ax4.transAxes)
                            plt.text(.2,1.1-2.6*nstascale,'Myy = {:0.3e}'.format(myy),transform=ax4.transAxes)
                            plt.text(.2,1.1-2.8*nstascale,'Myz = {:0.3e}'.format(myz),transform=ax4.transAxes)
                            plt.text(.2,1.1-3.0*nstascale,'Mzz = {:0.3e}'.format(mzz),transform=ax4.transAxes)
                            plt.axis('off')
                        else:
                            ax4 = fig.add_subplot(nsta,ncol,ncol)
                            plt.text(.2,.8,'Depth = ' + str(depth),transform=ax4.transAxes)
                            plt.text(.2,.8-.2*nstascale,'Strike = ' + str(round(st1)) + ' ; ' + str(round(st2)),transform=ax4.transAxes)
                            plt.text(.2,.8-.4*nstascale,'Rake = ' + str(round(rk1)) + ' ; ' + str(round(rk2)),transform=ax4.transAxes)
                            plt.text(.2,.8-.6*nstascale,'Dip = ' + str(round(dp1)) + ' ; ' + str(round(dp2)),transform=ax4.transAxes)
                            plt.text(.2,.8-.8*nstascale,'Mo = {:0.2e}'.format(Motot),transform=ax4.transAxes)
                            plt.text(.2,.8-1.0*nstascale,'Mw = ' + str(round(Mw,2)),transform=ax4.transAxes)
                            plt.text(.2,.8-1.2*nstascale,'Percent DC = ' + str(round(pdc,2)),transform=ax4.transAxes)
                            plt.text(.2,.8-1.4*nstascale,'Percent CLVD = ' + str(round(pclvd,2)),transform=ax4.transAxes)
                            plt.text(.2,.8-1.6*nstascale,'Percent ISO = ' + str(round(piso,2)),transform=ax4.transAxes)
                            plt.text(.2,.8-1.8*nstascale,'Var. Red. = ' + str(round(VRtot,2)),transform=ax4.transAxes)
                            plt.text(.2,.8-2.0*nstascale,'Mxx = {:0.3e}'.format(mxx),transform=ax4.transAxes)
                            plt.text(.2,.8-2.2*nstascale,'Mxy = {:0.3e}'.format(mxy),transform=ax4.transAxes)
                            plt.text(.2,.8-2.4*nstascale,'Mxz = {:0.3e}'.format(mxz),transform=ax4.transAxes)
                            plt.text(.2,.8-2.6*nstascale,'Myy = {:0.3e}'.format(myy),transform=ax4.transAxes)
                            plt.text(.2,.8-2.8*nstascale,'Myz = {:0.3e}'.format(myz),transform=ax4.transAxes)
                            plt.text(.2,.8-3.0*nstascale,'Mzz = {:0.3e}'.format(mzz),transform=ax4.transAxes)
                            plt.axis('off')





                    #plotting synthetic data
                    #it plots this late in the program because scaley doesn't work until the last plots
                    ax3.plot(ST,'r--',linewidth=.75, scaley=False)

                    ax6.plot(SR,'r--',linewidth=.75, scaley=False)

                    ax2.plot(SZ,'r--',linewidth=.75, scaley=False)

                    index=index+1


                    if truensta==index:
                        #creating beachball in large diagram
                        mt=np.array((mxx,myy,mzz,mxy,mxz,myz))
                        beach1 = beach(mt,xy=(.5,.5),width=.95,mopad_basis='NED',show_iso=True,facecolor='black') #Show Iso? Should it be true?
                        if nsta1==False:
                            if nsta==2:
                                ax5 = fig.add_subplot(3,ncol,3*ncol)
                            if nsta<=4 and nsta>2:
                                ax5 = fig.add_subplot(nsta,ncol,nsta*ncol)
                            if nsta>4:
                                ax5 = fig.add_subplot(4,ncol,4*ncol)
                        if nsta1==True:
                            ax5 = fig.add_subplot(nsta,ncol,nsta*ncol-2)
                        ax5.add_collection(beach1)
                        ax5.set_aspect("equal")
                        ax5.set_axis_off()

                        #plots the filenames on the beachball
                        if args.noplotazimuth==False:
                            for filenameaz in filenameazlist:
                                Az = (90 - (float((filenameaz.split())[0])))
                                if Az < 0:
                                    Az+=360
                                filename = (filenameaz.split())[1]
                                if Az>0 and Az<=90:
                                    ax5.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='left', va='bottom')
                                if Az>90 and Az<=180:
                                    ax5.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='right', va='bottom')
                                if Az>180 and Az<=270:
                                    ax5.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='right', va='top')
                                if Az>270 and Az<=360 or Az==0:
                                    ax5.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='left', va='top')

                    #plots and saves the current plot and makes a new page if there's more than 8 stations
                    if index==truensta or index==8:
                        if args.nosavefig==False:
                            pdf.savefig()
                        if args.nodisplaywindow==False:
                            fig.canvas.manager.set_window_title(mtinput)
                            plt.show()
                        plt.close()
                        if truensta>index:
                            fig=plt.figure(figsize=(9,6))
                            if (truensta-index)<8:
                                nsta=truensta-index
                                if nsta==1:
                                    nsta=2
                                    nsta1=True
                            else:
                                nsta=8





#                    """
#                    if index==nsta-1:
#                        # Plot beachball
#                        # In Obspy: M11, M22, M33, M12, M13, M23
#                        ax2=plt.subplot(nsta,ncol,4+index*ncol)
#                        ###plt.axis('equal')
#                        mt=[mxx,myy,mzz,mxy,mxz,myz]
#                        beach=beach(mt, xy=(0.5, 0.5), width=0.6)
#                        ax2.add_collection(beach)
#                    """

                #clear data list
                DT=[];DR=[];DZ=[]
                ST=[];SR=[];SZ=[]
                continue
            else:
                if mark==1:
                    depth=float((line.split())[0])
                    variance=float((line.split())[1])
                    VRtot=float((line.split())[2])
                    nsta=int((line.split())[3])
                    npages=math.ceil(nsta/8)
                    truensta=nsta
                    if nsta>8:
                        nsta=8
                    nsta1=False
                    if nsta==1:
                        nsta=2
                        nsta1=True
                    mark=0 #reset
                    continue
                elif mark==2:
                    mxx=float((line.split())[0])
                    mxy=float((line.split())[1])
                    mxz=float((line.split())[2])
                    myy=float((line.split())[3])
                    myz=float((line.split())[4])
                    mzz=float((line.split())[5])
                                                                                     #note if applied to tensors and passed to
                                                                                     #mopad the very large values leads to plotting errors can result
                    Mfull=np.array([[mxx,mxy,mxz],[mxy,myy,myz],[mxz,myz,mzz]])  #Construct Moment Tensor Matrix
                    L, V = la.eig(Mfull)

                    if L[0]==L[1] and L[0]==L[2]:
                        print('Pure Isotropic')                                         #deal with this perfect isotropic case
                        mxx=mxx+mxx*0.0001
                        Mfull[0,]=Mfull[0,]+Mfull[0,]*0.0001

                    Moiso=(mxx+myy+mzz)/3                                                #Esimate the Scalar Moment
                    Mdev=Mfull - np.identity(3)*Moiso                                    #Compute the Deviatoric Moment Tensor

                    w, v = la.eig(Mdev)                                                  #Calculate eigenvalues(w) and eigenvectors(v)

                    Motot=(abs(Moiso) + max(abs(w)))*moscale                            #Compute Bower and Hudson Total Moment and the
                    Mw=(np.log10(Motot)-16.1)/1.5                                      #Moment Magnitude

                    Moiso=Moiso*moscale                                                  #Now scale Moiso and Modev for plotting later
                    Modev=max(abs(w))*moscale                                      #Modev is maximum deviatoric eigenvalue in absolute sense                                                            #It is used to scale deviatoric tensor into DC and CLVD components


                    #Order the eigenvalues and eigenvectors
                    indx=np.argsort(abs(w))                                              #Sort by absolute value of w
                    m3=w[indx[2]]
                    m2=w[indx[1]]
                    m1=w[indx[0]]
                    eig3=v[:,indx[2]]
                    eig2=v[:,indx[1]]
                    eig1=v[:,indx[0]]

                    #Order eigenvalues for Tape & Tape Lune
                    indx=np.argsort(L)                                                  #Sort retaining sign
                    l1=L[indx[2]]
                    l2=L[indx[1]]
                    l3=L[indx[0]]

                    #Calculate Tape & Tape gamma and beta parameters testing for pure isotropic singularity
                    #These parameters, gamma, beta and delta are used later to plot the source-type in the Tape and Tape Lune perspective
                    if l1 == l2 and l1 == l3 and l1 > 0.:
                        gamma=0.
                        beta=0.
                        delta=90. - beta
                    elif l1 == l2 and l1 == l3 and l1 < 0.:
                        gamma=0.
                        beta=0.
                        delta=beta - 90.
                    else:
                        gamma=math.atan((-l1+2*l2-l3)/(np.sqrt(3)*(l1-l3)))*180/math.pi
                        beta=math.acos((l1+l2+l3)/(np.sqrt(3)*np.sqrt(L.dot(L))))*180/math.pi
                        delta=90. - beta

                    #Construct Dyadics
                    #Dyadics represent fundamental vector-dipole tensors from which double-couples, CLVDs, tensile-cracks, etc. are constructed
                    #See Jost and Herrman for details
                    a3=np.array((eig3, eig3, eig3)).transpose()
                    a2=np.array((eig2, eig2, eig2)).transpose()
                    a1=np.array((eig1, eig1, eig1)).transpose()
                    a3a3=a3*a3.transpose()
                    a2a2=a2*a2.transpose()
                    a1a1=a1*a1.transpose()

                    #Perform DC-CLVD Decomposition
                    F=-1*m1/m3
                    Mdc=m3*(1-2*F)*(a3a3-a2a2)                                  #Double-Couple Moment Tensor
                    Mclvd=m3*F*(2*a3a3-a2a2-a1a1)                               #CLVD Moment Tensor
                    Modc=abs(m3*(1-2*F))*moscale                                #Double-Couple Moment
                    Moclvd=abs(2*m3*F)*moscale                                  #CLVD Moment - to be consistent with Hudson decomp
                    kappa=Moiso/Motot                                           #Hudson Plot kappa
                    T=(2*m1)/abs(m3)                                            #Hudson Plot T
                    periso=abs(Moiso/Motot)
                    perdc=abs(Modc/Modev)
                    perclvd=abs(Moclvd/Modev)

                    #Determine Strike, Rake, Dip
                    if Modc != 0.:
                        w, v = la.eig(Mdc)
                        indx=np.argsort(w)                                              #Sort by absolute value of w
                        eig3=v[:,indx[2]]
                        eig2=v[:,indx[1]]
                        eig1=v[:,indx[0]]
                        nu1=(1/np.sqrt(2))*(eig3-eig1)   #fault normal vector
                        u1=(1/np.sqrt(2))*(eig1+eig3)    #slip vector
                        [strike1, rake1, dip1]=fpsol.fpsol(nu1,u1)
                        nu2=(1/np.sqrt(2))*(eig1+eig3)   #conjugate fault normal vector
                        u2=(1/np.sqrt(2))*(eig3-eig1)    #conjugate slip vector
                        [strike2, rake2, dip2]=fpsol.fpsol(nu2,u2)

                    mark=0
                    continue
                elif mark==3:
                    st1=float((line.split())[0])
                    rk1=float((line.split())[1])
                    dp1=float((line.split())[2])
                    st2=float((line.split())[3])
                    rk2=float((line.split())[4])
                    dp2=float((line.split())[5])
                    mark=0
                    continue
                elif mark==4:
                    pdc=float((line.split())[0])
                    pclvd=float((line.split())[1])
                    piso=float((line.split())[2])
                    mark=0 #reset
                    continue
                elif mark==5:
                    filename=line
                    #strips anything after the first '_' or '.' in the filename
                    terminator = filename.find('_')
                    if terminator==-1:
                        terminator = filename.find('.')
                    if terminator!=-1:
                        filename = filename[:terminator]
                    mark=0
                    continue
                elif mark==6:
                    dt=float((line.split())[0])
                    npts=int((line.split())[1])
                    dist=float((line.split())[2])
                    Az=float((line.split())[3])
                    filenameazlist.append(str(Az) + " " + filename)
                    Zcor=int((line.split())[4])
                    VR=float((line.split())[5])
                    mark=0
                    continue
                elif mark==7:
                    DT.append(float((line.split())[0]))
                    DR.append(float((line.split())[1]))
                    DZ.append(float((line.split())[2]))
                    ST.append(float((line.split())[3]))
                    SR.append(float((line.split())[4]))
                    SZ.append(float((line.split())[5]))

            # Locate yourself...
            if line.startswith('#depth'):
                mark=1;
            elif line.startswith('#mxx'):
                mark=2
            elif line.startswith('#st1'):
                mark=3
            elif line.startswith('#pdc'):
                mark=4
            elif line.startswith('#filename'):
                mark=5
            elif line.startswith('#dt'):
                mark=6
            elif line.startswith('#data'):
                mark=7

        f.close()


    if args.plotseparatebeachball==True:
        fig=plt.figure(figsize=(30,30))
        mt=np.array((mxx,myy,mzz,mxy,mxz,myz))
        #ax2=beachball(mt,facecolor='k',outfile='fm_beachball.pdf')
        beach1 = beach(mt,xy=(0.5,0.5),width=0.95,mopad_basis='NED',show_iso=True,facecolor='black')
        ax2 = fig.add_subplot(1,1,1)
        ax2.add_collection(beach1)
        ax2.set_aspect("equal")
        ax2.set_axis_off()
        if args.nosavefig==False:
            if args.noplotazimuth==False:
                for filenameaz in filenameazlist:
                    Az = 90 - float((filenameaz.split())[0])
                    if Az < 0:
                        Az+=360
                    filename = (filenameaz.split())[1]
                    if Az>0 and Az<=90:
                        ax2.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='left', va='bottom', fontsize=70)
                    if Az>90 and Az<=180:
                        ax2.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='right', va='bottom', fontsize=70)
                    if Az>180 and Az<=270:
                        ax2.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='right', va='top', fontsize=70)
                    if Az>270 and Az<=360 or Az==0:
                        ax2.text((.5+0.475*math.cos(np.deg2rad(Az))), (.5+0.475*math.sin(np.deg2rad(Az))), filename, ha='left', va='top', fontsize=70)
            fig.savefig('fm_beachball(' + mtinput + ').pdf')
        plt.close()



    if args.plottape==True:
        nssfilename=''
        nssplotflag=0

        Mfull=np.array([[mxx,mxy,mxz],[mxy,myy,myz],[mxz,myz,mzz]])  #Construct Moment Tensor Matrix
        L, V = la.eig(Mfull)

        if L[0]==L[1] and L[0]==L[2]:
            print('Pure Isotropic')                                         #deal with this perfect isotropic case
            mxx=mxx+mxx*0.0001
            Mfull[0,]=Mfull[0,]+Mfull[0,]*0.0001

        Moiso=(mxx+myy+mzz)/3                                                #Esimate the Scalar Moment

        Mdev=Mfull - np.identity(3)*Moiso                                    #Compute the Deviatoric Moment Tensor

        w, v = la.eig(Mdev)                                                  #Calculate eigenvalues(w) and eigenvectors(v)

        Motot=(abs(Moiso) + max(abs(w)))*moscale                             #Compute Bower and Hudson Total Moment and the
        Mw=(np.log10(Motot)-9.1)/1.5                                         #Moment Magnitude

        Moiso=Moiso*moscale                                                  #Now scale Moiso and Modev for plotting later
        Modev=max(abs(w))*moscale                                            #Modev is maximum deviatoric eigenvalue in absolute sense
                                                                             #It is used to scale deviatoric tensor into DC and CLVD components


        #Order the eigenvalues and eigenvectors
        indx=np.argsort(abs(w))                                              #Sort by absolute value of w
        m3=w[indx[2]]
        m2=w[indx[1]]
        m1=w[indx[0]]
        eig3=v[:,indx[2]]
        eig2=v[:,indx[1]]
        eig1=v[:,indx[0]]

        #Order eigenvalues for Tape & Tape Lune
        indx=np.argsort(L)                                                  #Sort retaining sign
        l1=L[indx[2]]
        l2=L[indx[1]]
        l3=L[indx[0]]

        #Calculate Tape & Tape gamma and beta parameters testing for pure isotropic singularity
        #These parameters, gamma, beta and delta are used later to plot the source-type in the Tape and Tape Lune perspective
        if l1 == l2 and l1 == l3 and l1 > 0.:
            gamma=0.
            beta=0.
            delta=90. - beta
        elif l1 == l2 and l1 == l3 and l1 < 0.:
            gamma=0.
            beta=0.
            delta=beta - 90.
        else:
            gamma=math.atan((-l1+2*l2-l3)/(np.sqrt(3)*(l1-l3)))*180/math.pi
            beta=math.acos((l1+l2+l3)/(np.sqrt(3)*np.sqrt(L.dot(L))))*180/math.pi
            delta=90. - beta

        #Construct Dyadics
        #Dyadics represent fundamental vector-dipole tensors from which double-couples, CLVDs, tensile-cracks, etc. are constructed
        #See Jost and Herrman for details
        a3=np.array((eig3, eig3, eig3)).transpose()
        a2=np.array((eig2, eig2, eig2)).transpose()
        a1=np.array((eig1, eig1, eig1)).transpose()
        a3a3=a3*a3.transpose()
        a2a2=a2*a2.transpose()
        a1a1=a1*a1.transpose()

        #Perform DC-CLVD Decomposition
        F=-1*m1/m3
        Mdc=m3*(1-2*F)*(a3a3-a2a2)                                  #Double-Couple Moment Tensor
        Mclvd=m3*F*(2*a3a3-a2a2-a1a1)                               #CLVD Moment Tensor
        Modc=abs(m3*(1-2*F))*moscale                                #Double-Couple Moment
        Moclvd=abs(2*m3*F)*moscale                                  #CLVD Moment - to be consistent with Hudson decomp
        kappa=Moiso/Motot                                           #Hudson Plot kappa
        T=(2*m1)/abs(m3)                                            #Hudson Plot T
        periso=abs(Moiso/Motot)
        perdc=abs(Modc/Modev)
        perclvd=abs(Moclvd/Modev)

        #Determine Strike, Rake, Dip
        if Modc != 0.:
            w, v = la.eig(Mdc)
            indx=np.argsort(w)                                              #Sort by absolute value of w
            eig3=v[:,indx[2]]
            eig2=v[:,indx[1]]
            eig1=v[:,indx[0]]
            nu1=(1/np.sqrt(2))*(eig3-eig1)   #fault normal vector
            u1=(1/np.sqrt(2))*(eig1+eig3)    #slip vector
            [strike1, rake1, dip1]=fpsol.fpsol(nu1,u1)
            nu2=(1/np.sqrt(2))*(eig1+eig3)   #conjugate fault normal vector
            u2=(1/np.sqrt(2))*(eig3-eig1)    #conjugate slip vector
            [strike2, rake2, dip2]=fpsol.fpsol(nu2,u2)

        #Construct Moment Tensor arrays for plotting
        fm=np.array((mxx,myy,mzz,mxy,mxz,myz))
        devm=np.array((Mdev[0,0], Mdev[1,1], Mdev[2,2], Mdev[0,1], Mdev[0,2], Mdev[1,2]))
        dcm=np.array((Mdc[0,0], Mdc[1,1], Mdc[2,2], Mdc[0,1], Mdc[0,2], Mdc[1,2]))
        clvdm=np.array((Mclvd[0,0], Mclvd[1,1], Mclvd[2,2], Mclvd[0,1], Mclvd[0,2], Mclvd[1,2]))

        #Compute Tape and Tape parameters for plotting the NSS
        #Read the NSS output for plotting solution space on Tape and Tape Lune
        if (nssplotflag == 1):
            data=pd.read_csv(nssfilename, sep='\s+', header=None)

            d=np.array(data)
            lam=np.array((d[:,0],d[:,1],d[:,2])).transpose()   #eigenvalues are column ordered each row is a individual tuple
            lam.sort(axis=1)   #sort eigenvalue rows lam1=d[:,2], lam2=d[:,1], lam3=d[:,0]
            vr=d[:,3]

            l1=lam[:,2]
            l2=lam[:,1]
            l3=lam[:,0]
            L=np.sqrt(l1**2 + l2**2 + l3**2)

        #Test for pure isotropic singularity and compute gamma, beta and delta
            n=len(l1)
            GAMMA=np.zeros(n)
            BETA=np.zeros(n)
            DELTA=np.zeros(n)
            for i in range(0,n,1):
                if l1[i] == l2[i] and l1[i] == l3[i] and l1[i] > 0.:
                    GAMMA[i]=0.
                    BETA[i]=0.
                    DELTA[i]=90. - BETA[i]
                elif l1[i] == l2[i] and l1[i] == l3[i] and l1[i] < 0.:
                    GAMMA[i]=0.
                    BETA[i]=0.
                    DELTA[i]=BETA[i] - 90.
                else:
                    GAMMA[i]=np.arctan((-l1[i]+2*l2[i]-l3[i])/(np.sqrt(3)*(l1[i]-l3[i])))*180/np.pi
                    BETA[i]=np.arccos((l1[i]+l2[i]+l3[i])/(np.sqrt(3)*L[i]))*180/np.pi
                    DELTA[i]=90. - BETA[i]
        #Plot Tape and Tape Lune
        #Initial code from 'Ajean' https://stackoverflow.com/questions/32209496/matplotlib-basemap-fundamental-lune

        import cartopy.crs as ccrs
        import matplotlib.path as mpath
        from scipy.interpolate import griddata

        # Mollweide projection
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(111, projection=ccrs.LambertAzimuthalEqualArea())  #This seems best
        ax.set_extent([-30, 30, -90, 90])
        xstep=30/5 #20% lines
        ystep=90/5 #20% lines
        xgrds=np.arange(-30.0, 31.0, xstep)
        ygrds=np.arange(-90.0, 91.0, ystep)
        ax.gridlines(xlocs=xgrds,ylocs=ygrds)

        # Here I define a matplotlib Path object to use as the boundary
        outlinex = np.concatenate([[-30],np.tile(-30,180), np.tile(30,180),[-30]])
        outliney = np.concatenate([[-90],np.arange(-90,90),np.arange(89,-91,-1),[-90]])
        outlinecodes = np.array([mpath.Path.MOVETO]+[mpath.Path.LINETO]*360+[mpath.Path.MOVETO])
        outlinepath = mpath.Path(np.column_stack([outlinex[::-1], outliney[::-1]]), outlinecodes[::-1])
        ax.set_boundary(outlinepath, transform=ccrs.Geodetic())



        #Fundamental Source-Types
        ax.plot(0, 90., 'ro', markersize=10, transform=ccrs.Geodetic())  #Explosion
        ax.text(30,87,'Explosion',fontsize=12,transform=ccrs.Geodetic())
        ax.plot(0, -90., 'ro', markersize=10, transform=ccrs.Geodetic())  #Implosion
        ax.text(70,-88,'Implosion',fontsize=12,transform=ccrs.Geodetic())
        ax.plot(0, 0, 'ro', markersize=10, transform=ccrs.Geodetic())   #Double-Couple
        ax.text(0,2,'DC',fontsize=12,transform=ccrs.Geodetic())
        ax.plot(30, 0, 'ro', markersize=10, transform=ccrs.Geodetic())  #Negative CLVD
        ax.text(31,0,'-CLVD',fontsize=12,transform=ccrs.Geodetic())
        ax.plot(-30, 0, 'ro', markersize=10, transform=ccrs.Geodetic()) #Positive CLVD
        ax.text(-39,0,'+CLVD',fontsize=12,transform=ccrs.Geodetic())
        LAM=np.array([3,1,1])
        x=math.atan((-LAM[0]+2*LAM[1]-LAM[2])/(np.sqrt(3)*(LAM[0]-LAM[2])))*180/math.pi
        y=math.acos((LAM[0]+LAM[1]+LAM[2])/(np.sqrt(3)*np.sqrt(LAM.dot(LAM))))*180/math.pi
        y=90. - y
        ax.plot(x, y, 'ro', markersize=10, transform=ccrs.Geodetic())  #Tensile Crack
        ax.text(x-15,y-2,'+Crack',fontsize=12,transform=ccrs.Geodetic())
        LAM=np.array([-1,-1, -3]) #note ordering is due to sign considered ordering
        x=math.atan((-LAM[0]+2*LAM[1]-LAM[2])/(np.sqrt(3)*(LAM[0]-LAM[2])))*180/math.pi
        y=math.acos((LAM[0]+LAM[1]+LAM[2])/(np.sqrt(3)*np.sqrt(LAM.dot(LAM))))*180/math.pi
        y=90. - y
        ax.plot(x, y, 'ro', markersize=10, transform=ccrs.Geodetic())  #Closing Crack
        ax.text(x+3,y-1,'-Crack',fontsize=12,transform=ccrs.Geodetic())
        LAM=np.array([1,0,0])
        x=math.atan((-LAM[0]+2*LAM[1]-LAM[2])/(np.sqrt(3)*(LAM[0]-LAM[2])))*180/math.pi
        y=math.acos((LAM[0]+LAM[1]+LAM[2])/(np.sqrt(3)*np.sqrt(LAM.dot(LAM))))*180/math.pi
        y=90. - y
        ax.plot(x, y, 'ro', markersize=10, transform=ccrs.Geodetic())  #LVD
        ax.text(x-10,y-2,'+LVD',fontsize=12,transform=ccrs.Geodetic())
        LAM=np.array([0,0,-1])
        x=math.atan((-LAM[0]+2*LAM[1]-LAM[2])/(np.sqrt(3)*(LAM[0]-LAM[2])))*180/math.pi
        y=math.acos((LAM[0]+LAM[1]+LAM[2])/(np.sqrt(3)*np.sqrt(LAM.dot(LAM))))*180/math.pi
        y=90. - y
        ax.plot(x, y, 'ro', markersize=10, transform=ccrs.Geodetic())  #LVD
        ax.text(x+3,y-0,'-LVD',fontsize=12,transform=ccrs.Geodetic())



        # Plot Source-Type Solution Space
        # Comment out if not available
        if (nssplotflag == 1):
            c = plt.cm.plasma(np.arange(0.,100.,10.)/100)
            x=np.arange(-30.,31,5)  #The third argument, the step controls smoothing
            y=np.arange(-90,90,5)
            X, Y= np.meshgrid(x, y)
            idx=np.nonzero(vr >= 10.)
            Z = griddata((GAMMA[idx],DELTA[idx]),vr[idx],(X,Y), method='cubic')
            cb=ax.contourf(X, Y, Z, 20, transform=ccrs.PlateCarree(),cmap='Blues')

        #Plot MT Solution
        ax.plot(gamma, delta, 'ks', markersize=14, transform=ccrs.Geodetic())

        # Add colorbar, make sure to specify tick locations to match desired ticklabels
        #ax.set_title('Source-Type Lune')
        if (nssplotflag == 1):
            position=fig.add_axes([0.70,0.3,0.025,0.4])  ## the parameters are the specified position you set
            cbar=plt.colorbar(cb, cax=position, orientation='vertical',ticks=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100], spacing='uniform',shrink=0.5)
            cbar.set_label('Variance Reduction (%)', rotation=90, size=14)

        fig.savefig('fm_tape(' + mtinput +').pdf')


    #Beginning of third section of code from Moment Tensor Decomposition Tool

    #Construct Moment Tensor arrays for plotting
    mt=np.array((mxx,myy,mzz,mxy,mxz,myz))
    devm=np.array((Mdev[0,0], Mdev[1,1], Mdev[2,2], Mdev[0,1], Mdev[0,2], Mdev[1,2]))
    dcm=np.array((Mdc[0,0], Mdc[1,1], Mdc[2,2], Mdc[0,1], Mdc[0,2], Mdc[1,2]))
    clvdm=np.array((Mclvd[0,0], Mclvd[1,1], Mclvd[2,2], Mclvd[0,1], Mclvd[0,2], Mclvd[1,2]))

    fig=plt.figure(figsize=(8,8))
    threshold=0.;                    #initialize threshold
    if Moiso != 0.0:
        beach1 = beach(mt,xy=(0.5,0.5),width=0.95,mopad_basis='NED',show_iso=True,facecolor='black')
        ax2 = fig.add_subplot(2,2,1)
        ax2.add_collection(beach1)
        ax2.set_aspect("equal")
        ax2.set_axis_off()
        buf="Full MT {0:.2e}".format(Motot)
        ax2.set(title=buf)
        threshold=Moiso*0.00001  #Set Modev threshold to a small value of Mosio if there is a Moiso

    # Second one plots deviatoric mt
    if Modev != 0.0 and Modev/(Motot) > 0.001:   #plot only significant deviatoric parts
        beach1 = beach(devm,xy=(0.5,0.5),width=0.95*Modev/Motot,mopad_basis='NED',facecolor='black')
        ax3 = fig.add_subplot(2,2,2)
        ax3.add_collection(beach1)
        ax3.set_aspect("equal")
        ax3.set_axis_off()
        buf="Dev MT {0:.2e}".format(Modev)
        ax3.set(title=buf)

    # Third one plots dc
    if Modc != 0.0 and Modc/Motot > 0.001:     #plot only significant double-couple parts
        beach1 = beach(dcm,xy=(0.5,0.5),width=0.95*Modc/Modev,mopad_basis='NED',facecolor='black')
        ax3 = fig.add_subplot(2,2,3)
        ax3.add_collection(beach1)
        ax3.set_aspect("equal")
        ax3.set_axis_off()
        buf="DC MT {0:.2e}".format(Modc)
        ax3.set(title=buf)

    # Forth one plots dc
    if Moclvd != 0.0 and Moclvd/Motot > 0.001:  #plot only signicant clvd parts
        beach1 = beach(clvdm,xy=(0.5,0.5),width=0.95*Moclvd/Modev,mopad_basis='NED',facecolor='black')
        ax3 = fig.add_subplot(2,2,4)
        ax3.add_collection(beach1)
        ax3.set_aspect("equal")
        ax3.set_axis_off()
        buf="CLVD MT {0:.2e}".format(Moclvd)
        ax3.set(title=buf)

    if args.nosavefig==False:
        fig.savefig('fm_beachball_components(' + mtinput + ').pdf')
    #plt.show()
    plt.close()
