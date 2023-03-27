import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout
import numpy as np
import pandas as pd
import os
from matplotlib import cm
import matplotlib.dates as dts
import xarray as xr
import sys

class NaisChecker(QWidget):
    def __init__(self,data_file,boundary_file):
        super().__init__()
        
        self.window = QWidget()
        
        self.data_file = data_file
        self.boundary_file = boundary_file
        
        figWidget = pg.GraphicsLayoutWidget()
        layout = QVBoxLayout()
        
        self.negIonPlot = figWidget.addPlot(row=0,col=0,title="Negative ions")
        self.posIonPlot = figWidget.addPlot(row=1,col=0,title="Positive ions")
        self.ionFlagPlot = figWidget.addPlot(row=2,col=0,title="Ion flags")
        self.negParPlot = figWidget.addPlot(row=3,col=0,title="Negative particles")
        self.posParPlot = figWidget.addPlot(row=4,col=0,title="Positive particles")        
        self.parFlagPlot = figWidget.addPlot(row=5,col=0,title="Particle flags")
        
        boundButton = QPushButton("Save boundaries")
        boundButton.setFont(QtGui.QFont("sans",20))
        boundButton.setStyleSheet("background-color : lightgrey")
        boundButton.clicked.connect(lambda: self.saveBoundaries())
    
        layout.addWidget(figWidget)
        layout.addWidget(boundButton)
        self.window.setLayout(layout)
                
        self.window.show()
        
        # Create plot items from nais data
        self.createNaisImg()
    
        self.negIonPlot.addItem(self.negIonImg)
        self.posIonPlot.addItem(self.posIonImg)
        self.negParPlot.addItem(self.negParImg)
        self.posParPlot.addItem(self.posParImg)
        
        self.ionFlagPlot.addItem(self.ionScatter)
        self.parFlagPlot.addItem(self.parScatter)
        
        ionAx = self.ionFlagPlot.getAxis("left")
        parAx = self.parFlagPlot.getAxis("left")
        
        ionAx.setTicks([self.negIonFlagTicks + self.posIonFlagTicks])
        parAx.setTicks([self.negParFlagTicks + self.posParFlagTicks])
        
        self.ionFlagPlot.showGrid(y=True)
        self.parFlagPlot.showGrid(y=True)
        
        # Link the axes for synchronized panning/zooming
        self.posIonPlot.setXLink(self.negIonPlot)
        self.negParPlot.setXLink(self.negIonPlot)
        self.posParPlot.setXLink(self.negIonPlot)
        self.ionFlagPlot.setXLink(self.negIonPlot)
        self.parFlagPlot.setXLink(self.negIonPlot)

        self.posIonPlot.setYLink(self.negIonPlot)
        self.negParPlot.setYLink(self.negIonPlot)
        self.posParPlot.setYLink(self.negIonPlot)
     
        # Create adjustable colorbar for each plot
        cmap = pg.colormap.getFromMatplotlib('turbo')
    
        negIonColorbar = pg.ColorBarItem(colorMap=cmap)
        posIonColorbar = pg.ColorBarItem(colorMap=cmap)
        negParColorbar = pg.ColorBarItem(colorMap=cmap)
        posParColorbar = pg.ColorBarItem(colorMap=cmap)
        
        negIonColorbar.setImageItem(self.negIonImg,insert_in=self.negIonPlot)
        posIonColorbar.setImageItem(self.posIonImg,insert_in=self.posIonPlot)
        negParColorbar.setImageItem(self.negParImg,insert_in=self.negParPlot)
        posParColorbar.setImageItem(self.posParImg,insert_in=self.posParPlot)
    
        negIonColorbar.setLevels([1, 4])
        posIonColorbar.setLevels([1, 4])
        negParColorbar.setLevels([1, 5])
        posParColorbar.setLevels([1, 5])
        
        # Create ROI with double click
        figWidget.scene().sigMouseClicked.connect(self.onClick)
        
        
    def onClick(self,event):      
        if ((event.button()==1) & event.double()):
            
            p_neg_ion = self.negIonPlot.getViewBox().mapSceneToView(event.scenePos())
            p_pos_ion = self.posIonPlot.getViewBox().mapSceneToView(event.scenePos())
            p_neg_par = self.negParPlot.getViewBox().mapSceneToView(event.scenePos())
            p_pos_par = self.posParPlot.getViewBox().mapSceneToView(event.scenePos())
            negIonRect = self.negIonPlot.getViewBox().viewRect()
            posIonRect = self.posIonPlot.getViewBox().viewRect()
            negParRect = self.negParPlot.getViewBox().viewRect()
            posParRect = self.posParPlot.getViewBox().viewRect()
            
            if negIonRect.contains(p_neg_ion):
                self.addNegIonRoi((p_neg_ion.x(),p_neg_ion.y()),(0,0))            
            elif posIonRect.contains(p_pos_ion):
                self.addPosIonRoi((p_pos_ion.x(),p_pos_ion.y()),(0,0)) 
            elif negParRect.contains(p_neg_par):
                self.addNegParRoi((p_neg_par.x(),p_neg_par.y()),(0,0)) 
            elif negParRect.contains(p_pos_par):
                self.addPosParRoi((p_pos_par.x(),p_pos_par.y()),(0,0)) 
            else:
                pass
            
    def createNaisImg(self):
        ds = xr.open_dataset(self.data_file)
        
        # Data
        negIonData = np.log10(ds.neg_ions.where(ds.neg_ions>0).values)
        posIonData = np.log10(ds.pos_ions.where(ds.pos_ions>0).values)
        negParData = np.log10(ds.neg_particles.where(ds.neg_particles>0).values)
        posParData = np.log10(ds.pos_particles.where(ds.pos_particles>0).values)
        
        x = dts.date2num(ds.time.values.min())
        w = dts.date2num(ds.time.values.max())-x
        y = np.log10(ds.diameter).min()
        h = np.log10(ds.diameter).max()-y
        boundaries = (x,y,w,h)
        
        # Flags        
        negIonFlags = ds.neg_ion_flags.where(ds.neg_ion_flags==1,np.nan).dropna(dim="flag", how="all")
        negIonTicks = np.arange(len(negIonFlags.flag))
        negIonFlags = negIonFlags*negIonTicks
        self.negIonFlagTicks = list(zip(negIonTicks,negIonFlags.flag.values))
        
        posIonFlags = ds.pos_ion_flags.where(ds.pos_ion_flags==1,np.nan).dropna(dim="flag", how="all")
        posIonTicks = np.arange(len(negIonFlags.flag),len(negIonFlags.flag) + len(posIonFlags.flag))
        posIonFlags = posIonFlags*posIonTicks
        self.posIonFlagTicks = list(zip(posIonTicks,posIonFlags.flag.values))
        
        negParFlags = ds.neg_particle_flags.where(ds.neg_particle_flags==1,np.nan).dropna(dim="flag", how="all")
        negParTicks = np.arange(len(negParFlags.flag))
        negParFlags = negParFlags*negParTicks
        self.negParFlagTicks = list(zip(negParTicks,negParFlags.flag.values))

        posParFlags = ds.pos_particle_flags.where(ds.pos_particle_flags==1,np.nan).dropna(dim="flag", how="all")
        posParTicks = np.arange(len(negParFlags.flag),len(negParFlags.flag) + len(posParFlags.flag))
        posParFlags = posParFlags*posParTicks
        self.posParFlagTicks = list(zip(posParTicks,posParFlags.flag.values))
        
        ds.close()
     
        # Data
        self.negIonImg = pg.ImageItem(negIonData,rect=boundaries)
        self.posIonImg = pg.ImageItem(posIonData,rect=boundaries)
        self.negParImg = pg.ImageItem(negParData,rect=boundaries)
        self.posParImg = pg.ImageItem(posParData,rect=boundaries)
        
        # Flags
        self.ionScatter = pg.ScatterPlotItem()
        for i,flag in self.negIonFlagTicks:
            data = negIonFlags.sel(flag=flag)
            self.ionScatter.addPoints(dts.date2num(data.time.values),data.values,pen=None, symbol='o', symbolSize=5, brush="b")
        
        for i,flag in self.posIonFlagTicks:
            data = posIonFlags.sel(flag=flag)
            self.ionScatter.addPoints(dts.date2num(data.time.values),data.values,pen=None, symbol='o', symbolSize=5, brush="r")
        
        self.parScatter = pg.ScatterPlotItem()
        for i,flag in self.negParFlagTicks:
            data = negParFlags.sel(flag=flag)
            self.parScatter.addPoints(dts.date2num(data.time.values),data.values,pen=None, symbol='o', symbolSize=5, brush="b")
        
        for i,flag in self.posParFlagTicks:
            data = posParFlags.sel(flag=flag)
            self.parScatter.addPoints(dts.date2num(data.time.values),data.values,pen=None, symbol='o', symbolSize=5, brush="r")
        
    def addNegIonRoi(self,origin,size):
        roi = pg.RectROI(origin,size,removable=True)
        roi.sigRemoveRequested.connect(self.removeNegIonRoi)
        self.negIonPlot.addItem(roi)        
        
    def removeNegIonRoi(self,event):
        self.negIonPlot.removeItem(event)   
        
    def addPosIonRoi(self,origin,size):
        roi = pg.RectROI(origin,size,removable=True)
        roi.sigRemoveRequested.connect(self.removePosIonRoi)
        self.posIonPlot.addItem(roi)        
        
    def removePosIonRoi(self,event):
        self.posIonPlot.removeItem(event)   
        
    def addNegParRoi(self,origin,size):
        roi = pg.RectROI(origin,size,removable=True)
        roi.sigRemoveRequested.connect(self.removeNegParRoi)
        self.negParPlot.addItem(roi)        
        
    def removeNegParRoi(self,event):
        self.negParPlot.removeItem(event)   
        
    def addPosParRoi(self,origin,size):
        roi = pg.RectROI(origin,size,removable=True)
        roi.sigRemoveRequested.connect(self.removePosParRoi)
        self.posParPlot.addItem(roi)        
        
    def removePosParRoi(self,event):
        self.posParPlot.removeItem(event)   

    def saveBoundaries(self):                        
        plots=[self.negIonPlot,self.posIonPlot,self.negParPlot,self.posParPlot]        
        ds=xr.Dataset()
        
        for idx, plot in enumerate(plots):        
            roiItems = [item for item in plot.allChildItems() if type(item)==pg.graphicsItems.ROI.RectROI]        
            df = pd.DataFrame(columns = ["time_left","time_right","diam_bottom","diam_top"])        
            i=0
            
            for roi in roiItems:
                x_upper_left,y_upper_left = roi.pos()        
                width,height = roi.size()                
                time_left = x_upper_left
                time_right = x_upper_left+width
                diam_bottom = 10**(y_upper_left-height)
                diam_top = 10**y_upper_left                
                df.loc[i] = [time_left,time_right,diam_bottom,diam_top]                
                i+=1
                
            if idx==0:
                ds = ds.assign_coords(neg_ion_roi_id=df.index)
                ds = ds.assign(neg_ion_time_left=("neg_ion_roi_id",df["time_left"]))
                ds = ds.assign(neg_ion_time_right=("neg_ion_roi_id",df["time_right"]))
                ds = ds.assign(neg_ion_diam_left=("neg_ion_roi_id",df["diam_bottom"]))
                ds = ds.assign(neg_ion_diam_right=("neg_ion_roi_id",df["diam_top"]))                
            elif idx==1:
                ds = ds.assign_coords(pos_ion_roi_id=df.index)
                ds = ds.assign(pos_ion_time_left=("pos_ion_roi_id",df["time_left"]))
                ds = ds.assign(pos_ion_time_right=("pos_ion_roi_id",df["time_right"]))
                ds = ds.assign(pos_ion_diam_left=("pos_ion_roi_id",df["diam_bottom"]))
                ds = ds.assign(pos_ion_diam_right=("pos_ion_roi_id",df["diam_top"]))    
            elif idx==2:
                ds = ds.assign_coords(neg_par_roi_id=df.index)
                ds = ds.assign(neg_par_time_left=("neg_par_roi_id",df["time_left"]))
                ds = ds.assign(neg_par_time_right=("neg_par_roi_id",df["time_right"]))
                ds = ds.assign(neg_par_diam_left=("neg_par_roi_id",df["diam_bottom"]))
                ds = ds.assign(neg_par_diam_right=("neg_par_roi_id",df["diam_top"]))    
            else:
                ds=ds.assign_coords(pos_par_roi_id=df.index)
                ds = ds.assign(pos_par_time_left=("pos_par_roi_id",df["time_left"]))
                ds = ds.assign(pos_par_time_right=("pos_par_roi_id",df["time_right"]))
                ds = ds.assign(pos_par_diam_left=("pos_par_roi_id",df["diam_bottom"]))
                ds = ds.assign(pos_par_diam_right=("pos_par_roi_id",df["diam_top"]))

        ds=ds.assign_attrs(
            {"roi_ids": "Number the ROIs for each polarity/mode",
             "diams":"Particle diameter in meters",
             "times":"Days since 1970-01-01 UTC (matplotlib format)"}   
        )

        ds.to_netcdf(self.boundary_file)

        print("Saved:", self.boundary_file)


def startNaisChecker(dataset_path,bounding_boxes_path):
    """
    Manually check a NAIS dataset and draw bounding boxes around bad data
    
    Parameters
    ----------
    
    data_file : str
        Name of NAIS netcdf data file including path
    boundary_file : str
        Name of file where to save the coordinates 
        of bad data bounding boxes.
    """
    app = QApplication([])
    NaisChecker(dataset_path,bounding_boxes_path)
    app.exec_()
