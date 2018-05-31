# Written by Duncan Forgan, 17/07/12
# Module contains several functions to streamline plotting of many columns from a 2D array

import matplotlib.pyplot as plt
import numpy as np

def multigraph(data,nplots,xlabel, ylabels,setylog,ymin,ymax,outputstring):
    '''Plots many graphs from a 2D array'''
    # Python Function takes in several inputs

    # data - 2D array of data (rows,columns)
    # nplots - number of plots (from column 1 to nplots)
    # xlabel - string for the x axis
    # ylabels - array of strings for the y axis
    # setylog - boolean deciding whether y is log or not
    # ymin, ymax - arrays describing axis limits
    # outputstring - string determining the filenames to be outputted

    print 'Calling multigraph function'
    
    # Begin plotting loop

    for i in range(1,nplots):
        print 'Plotting ',outputstring[i]
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.set_xlabel(ylabels[0])    
        ax.set_ylabel(ylabels[i])
    
        if setylog[0]:
            ax.set_xscale('log') 
    
        if setylog and sum(data[:,i])!=0.0:
            ax.set_yscale('log')
    
        if ymin[i]!=ymax[i]:
            print 'Setting limits ',ymin[i],ymax[i]
            ax.set_ylim(ymin[i],ymax[i])
    
        ax.plot(data[:,0],data[:,i])
        outputfile = outputstring[i]+'.png'
        plt.savefig(outputfile, format='png')

    # Finish profile plots
    print 'Multigraph call complete'
    
def multigraph_legend(data,nplots,xlabel, ylabels,setylog,ymin,ymax,outputstring,legendstring):
    '''Plots multiple graphs from 2D arrays with a specified legend '''
    
    # Python Function takes in several inputs

    # data - 2D array of data (rows,columns)
    # nplots - number of plots (from column 1 to nplots)
    # xlabel - string for the x axis
    # ylabels - array of strings for the y axis
    # setylog - boolean deciding whether y is log or not
    # ymin, ymax - arrays describing axis limits
    # outputstring - string determining the filenames to be outputted
    # legendstring - string determining the legends for each plot

    print 'Calling multigraph_legend function'
    
    # Begin plotting loop

    for i in range(1,nplots):
        print 'Plotting ',outputstring[i]
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.set_xlabel(ylabels[0])
        ax.set_ylabel(ylabels[i])
        
        if setylog[0]:
            ax.set_xscale('log')        
    
        if setylog[i] and sum(data[:,i])!=0.0:
            ax.set_yscale('log')
    
        if ymin[i]!=ymax[i]:
            print 'Setting limits ',ymin[i],ymax[i]
            ax.set_ylim(ymin[i],ymax[i])
    
        ax.plot(data[:,0],data[:,i], label=legendstring[i])
        ax.legend()
        outputfile = outputstring[i]+'.png'
        plt.savefig(outputfile, format='png')

    # Finish profile plots
    print 'Multigraph_legend call complete'
    
def multigraph_legend_points(data,nplots,xlabel,ylabels,setylog,ymin,ymax,outputstring,legendstring, xpoints,ypoints,sizepoints,pointcolours):
    '''Plots multiple graphs from 2D arrays with a specified legend with extra points added (circle patch)'''
    
    # Python Function takes in several inputs

    # data - 2D array of data (rows,columns)
    # nplots - number of plots (from column 1 to nplots)
    # xlabel - string for the x axis
    # ylabels - array of strings for the y axis
    # setylog - boolean deciding whether y is log or not
    # ymin, ymax - arrays describing axis limits
    # outputstring - string determining the filenames to be outputted
    # legendstring - string determining the legends for each plot
    # xpoints - x locations of circles
    # ypoints - y locations of circles
    # sizepoints - sizes of circles 

    print 'Calling multigraph_legend_with_points function'
    
    # Begin plotting loop

    
    for i in range(1,nplots):
        print 'Plotting ',outputstring[i]
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.set_xlabel(ylabels[0])
        ax.set_ylabel(ylabels[i])
        
        if setylog[0]:
            ax.set_xscale('log')        
    
        if setylog[i] and sum(data[:,i])!=0.0:
            ax.set_yscale('log')
    
        if ymin[i]!=ymax[i]:
            print 'Setting limits ',ymin[i],ymax[i]
            ax.set_ylim(ymin[i],ymax[i])
                
        ax.plot(data[:,0],data[:,i], label=legendstring[i])
        
        yplot = np.zeros(len(xpoints))
        yplot[:] = ypoints[:,i]            
        ax.scatter(xpoints,yplot,s=100,facecolor=pointcolours)
                            
        ax.legend()
        outputfile = outputstring[i]+'.png'
        plt.savefig(outputfile, format='png')

    # Finish profile plots
    print 'Multigraph_legend_with_points call complete'
