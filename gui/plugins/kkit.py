import sys
import os
import math
import re
from PyQt4 import QtGui, QtCore, Qt
#import pygraphviz as pgv
import networkx as nx
sys.path.insert(0, '/home/harsha/async/gui')
import numpy as np
import config
import pickle 
from default import *
from moose import *
sys.path.append('plugins')
from mplugin import *
from kkitUtil import *
from kkitQGraphics import PoolItem, ReacItem,EnzItem,CplxItem,ComptItem
from kkitViewcontrol import *
from kkitCalcArrow import *
from kkitOrdinateUtil import *
import posixpath
from mtoolbutton import MToolButton

class KkitPlugin(MoosePlugin):
    """Default plugin for MOOSE GUI"""
    def __init__(self, *args):
        #print args
        MoosePlugin.__init__(self, *args)
        self.view = None

    def getPreviousPlugin(self):
        return None

    def getNextPlugin(self):
        return None

    def getAdjacentPlugins(self):
        return []

    def getViews(self):
        return self._views

    def getCurrentView(self):
        return self.currentView

    def getEditorView(self):
        if not hasattr(self, 'editorView'):
            self.editorView = KkitEditorView(self)
            self.editorView.getCentralWidget().editObject.connect(self.mainWindow.objectEditSlot)
            self.currentView = self.editorView
        return self.editorView

    def getRunView(self):
        if self.view is not None: return self.view
        self.view = RunView(self)
        graphView = self.view._centralWidget
        graphView.setDataRoot(self.modelRoot)
        schedulingDockWidget = self.view.getSchedulingDockWidget().widget()
        self._kkitWidget = self.view.plugin.getEditorView().getCentralWidget()
        self.runView = KkitRunView(self)
        self.currentRunView = self.runView.getCentralWidget()
        schedulingDockWidget.runner.update.connect(self.currentRunView.changeBgSize)
        schedulingDockWidget.runner.resetAndRun.connect(self.currentRunView.resetColor)
        graphView.layout().addWidget(self.currentRunView,0,0,2,1)
        return self.view

class KkitRunView(MooseEditorView):

    def __init__(self, plugin):
        MooseEditorView.__init__(self, plugin)
    '''
    def getToolPanes(self):
        return super(KkitRunView, self).getToolPanes()

    def getLibraryPane(self):
        return super(KkitRunView, self).getLibraryPane()

    def getOperationsWidget(self):
        return super(KkitRunView, self).getOperationsPane()

    def getToolBars(self):
        return self._toolBars
    '''
    def getCentralWidget(self):
        if self._centralWidget is None:
            self._centralWidget = kineticRunWidget()
            self._centralWidget.setModelRoot(self.plugin.modelRoot)
        return self._centralWidget


class KkitEditorView(MooseEditorView):
    def __init__(self, plugin):
        MooseEditorView.__init__(self, plugin)
        ''' EditorView and if kkit model is loaded then save model in XML is allowed '''
      
        #self.insertMenu = QtGui.QMenu('&Insert')
        self.fileinsertMenu = QtGui.QMenu('&File')
        if not hasattr(self,'SaveModelAction'):
            self.saveModelAction = QtGui.QAction('SaveToSBMLFile', self)
            self.saveModelAction.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+S", None, QtGui.QApplication.UnicodeUTF8))
            self.connect(self.saveModelAction, QtCore.SIGNAL('triggered()'), self.SaveModelDialogSlot)
            self.fileinsertMenu.addAction(self.saveModelAction)
        self._menus.append(self.fileinsertMenu)

 
    def SaveModelDialogSlot(self):
        type_sbml = 'SBML'
        filters = {'SBML(*.xml)': type_sbml}
        filename,filter_ = QtGui.QFileDialog.getSaveFileNameAndFilter(None,'Save File','',';;'.join(filters))
        extension = ""
        if str(filename).rfind('.') != -1:
            filename = filename[:str(filename).rfind('.')]
            if str(filter_).rfind('.') != -1:
                extension = filter_[str(filter_).rfind('.'):len(filter_)-1]
        if filename:
            filename = filename+extension
            if filters[str(filter_)] == 'SBML':
                moose.writeSBML(str(filename),self.plugin.modelRoot)
    '''
    def getToolPanes(self):
        return super(KkitEditorView, self).getToolPanes()

    def getLibraryPane(self):
        return super(KkitEditorView, self).getLibraryPane()

    def getOperationsWidget(self):
        return super(KkitEditorView, self).getOperationsPane()

    def getToolBars(self):
        return self._toolBars
    '''
    def getCentralWidget(self):
        if self._centralWidget is None:
            self._centralWidget = kineticEditorWidget()
            self._centralWidget.setModelRoot(self.plugin.modelRoot)
        return self._centralWidget


class  KineticsWidget(EditorWidgetBase):
    def __init__(self, *args):
        EditorWidgetBase.__init__(self, *args)
        self.setAcceptDrops(True)
        self.border = 10        
        self.sceneContainer = QtGui.QGraphicsScene(self)
        self.sceneContainer.setSceneRect(self.sceneContainer.itemsBoundingRect())
        self.sceneContainer.setBackgroundBrush(QtGui.QColor(230,220,219,120))
        #self.insertMenu = QtGui.QMenu('&Insert')
        #self._menus.append(self.insertMenu)
        #self.insertMapper = QtCore.QSignalMapper(self)
        classlist = ['CubeMesh','CylMesh','Pool','FuncPool','SumFunc','Reac','Enz','MMenz','StimulusTable','Table']
        #insertMapper, actions = self.getInsertActions(classlist)

        #for action in actions:
        #    self.insertMenu.addAction(action)
        #print "viewType",self.viewType
        # pickled the color map file """
        colormap_file = open(os.path.join(config.settings[config.KEY_COLORMAP_DIR], 'rainbow2.pkl'),'rb')
        self.colorMap = pickle.load(colormap_file)
        colormap_file.close()   
        
        if self.modelRoot == '/':
            self.m = wildcardFind('/##[ISA=ChemCompt]')
        else:
            self.m = wildcardFind(self.modelRoot+'/##[ISA=ChemCompt]')
        if self.m:
            self.xmin = 0.0
            self.xmax = 1.0
            self.ymin = 0.0
            self.ymax = 1.0
            self.autoCordinatepos = {}
            
            # Compartment and its members are setup
            self.meshEntry,self.xmin,self.xmax,self.ymin,self.ymax,self.noPositionInfo = setupMeshObj(self.modelRoot)
            self.autocoordinates = False
            # srcdesConnection dictonary will have connection information between src and des """
            self.srcdesConnection = {}
            setupItem(self.modelRoot,self.srcdesConnection)

            if self.noPositionInfo:
                self.autocoordinates = True
                QtGui.QMessageBox.warning(self, 
                                          'No coordinates found for the model', 
                                          '\n Automatic layouting will be done')
            #raise Exception('Unsupported kkit version')
                self.xmin,self.xmax,self.ymin,self.ymax,self.autoCordinatepos = autoCoordinates(self.meshEntry,self.srcdesConnection)
            # TODO: size will be dummy at this point, but I need the availiable size from the Gui
            self.size = QtCore.QSize(1024 ,768)
            #self.size = QtCore.QSize(300,400)
            
            # Scale factor to translate the x -y position to fit the Qt graphicalScene, scene width. """
            if self.xmax-self.xmin != 0:
                self.xratio = (self.size.width()-10)/(self.xmax-self.xmin)
            else: self.xratio = self.size.width()-10
            
            if self.ymax-self.ymin:
                self.yratio = (self.size.height()-10)/(self.ymax-self.ymin)
            else: self.yratio = (self.size.height()-10)

        #A map b/w moose compartment key with QGraphicsObject
            self.qGraCompt = {}
        
        #A map between mooseId of all the mooseObject (except compartment) with QGraphicsObject
            self.mooseId_GObj = {}
        
            self.border = 5
            self.arrowsize = 2
            self.iconScale = 1
            self.defaultComptsize = 5
            self.itemignoreZooming = False
            self.lineItem_dict = {}
            self.object2line = defaultdict(list)

    def makePoolItem(self, poolObj, qGraCompt):
        raise NotImplementedError('method must be reimplemented in subclass')

    def mooseObjOntoscene(self):
        #  All the compartments are put first on to the scene \
        #  Need to do: Check With upi if empty compartments exist 
        for cmpt in sorted(self.meshEntry.iterkeys()):
            self.createCompt(cmpt)
            self.qGraCompt[cmpt]
            #comptRef = self.qGraCompt[cmpt]
        
        #Enzymes of all the compartments are placed first, \
        #     so that when cplx (which is pool object) queries for its parent, it gets its \
        #     parent enz co-ordinates with respect to QGraphicsscene """
        
        for cmpt,memb in self.meshEntry.items():
            for enzObj in find_index(memb,'enzyme'):
                enzinfo = enzObj.path+'/info'
                if enzObj.className == 'Enz':
                    enzItem = EnzItem(enzObj,self.qGraCompt[cmpt])
                else:
                    enzItem = MMEnzItem(enzObj,self.qGraCompt[cmpt])
                self.setupDisplay(enzinfo,enzItem,"enzyme")
                self.setupSlot(enzObj,enzItem)
        
        for cmpt,memb in self.meshEntry.items():
            for poolObj in find_index(memb,'pool'):
                poolinfo = poolObj.path+'/info'
                ''' depending on Editor Widget or Run widget pool will be created a PoolItem or PoolItemCircle '''
                poolItem = self.makePoolItem(poolObj,self.qGraCompt[cmpt])
                self.setupDisplay(poolinfo,poolItem,"pool")
                self.setupSlot(poolObj,poolItem)
            
            for cplxObj in find_index(memb,'cplx'):
                cplxinfo = (cplxObj.parent).path+'/info'
                cplxItem = CplxItem(cplxObj,self.mooseId_GObj[element(cplxObj).parent.getId()])
                self.setupDisplay(cplxinfo,cplxItem,"cplx")
                self.setupSlot(cplxObj,cplxItem)
            
            for reaObj in find_index(memb,'reaction'):
                reainfo = reaObj.path+'/info'
                reaItem = ReacItem(reaObj,self.qGraCompt[cmpt])
                self.setupDisplay(reainfo,reaItem,"reaction")
                self.setupSlot(reaObj,reaItem)
                
            for tabObj in find_index(memb,'table'):
                tabinfo = tabObj.path+'/info'
                tabItem = TableItem(tabObj,self.qGraCompt[cmpt])
                self.setupDisplay(tabinfo,tabItem,"tab")
                self.setupSlot(tabObj,tabItem)
        
        # compartment's rectangle size is calculated depending on children 
        for k, v in self.qGraCompt.items():
            rectcompt = v.childrenBoundingRect()
            v.setRect(rectcompt.x()-10,rectcompt.y()-10,(rectcompt.width()+20),(rectcompt.height()+20))
            v.setPen(QtGui.QPen(Qt.QColor(66,66,66,100), 5, Qt.Qt.SolidLine, Qt.Qt.RoundCap, Qt.Qt.RoundJoin))
            v.cmptEmitter.connect(v.cmptEmitter,QtCore.SIGNAL("qgtextPositionChange(PyQt_PyObject)"),self.positionChange)
            v.cmptEmitter.connect(v.cmptEmitter,QtCore.SIGNAL("qgtextItemSelectedChange(PyQt_PyObject)"),self.objectEditSlot)
        
    def createCompt(self,key):
        self.new_Compt = ComptItem(self,0,0,0,0,key)
        self.qGraCompt[key] = self.new_Compt
        self.new_Compt.setRect(10,10,10,10)
        self.sceneContainer.addItem(self.new_Compt)
        
    def setupDisplay(self,info,graphicalObj,objClass):
        xpos,ypos = self.positioninfo(info)
        # For Reaction and Complex object I have skipped the process to get the facecolor and background color as \
        #    we are not using these colors for displaying the object so just passing dummy color white 

        if( (objClass == "reaction" ) or (objClass == "cplx")):
            textcolor,bgcolor = "white","white"
        elif objClass == "tab":
            textcolor,bgcolor = getColor(info,self.colorMap)
        else:
            textcolor,bgcolor = getColor(info,self.colorMap)
        graphicalObj.setDisplayProperties(xpos,ypos,textcolor,bgcolor)
    
    def positioninfo(self,iteminfo):
        if self.noPositionInfo:
            
            try:
                # kkit does exist item's/info which up querying for parent.path gives the information of item's parent 
                x,y = self.autoCordinatepos[(element(iteminfo).parent).path]
            except:
                # But in Cspace reader doesn't create item's/info, up on querying gives me the error which need to change\
                # in ReadCspace.cpp, at present i am taking care b'cos i don't want to pass just the item where I need to check\
                # type of the object (rea,pool,enz,cplx,tab) which I have already done. 
                parent, child = posixpath.split(iteminfo)
                x,y = self.autoCordinatepos[parent]
            ypos = (y-self.ymin)*self.yratio
        else:
            x = float(element(iteminfo).getField('x'))
            y = float(element(iteminfo).getField('y'))
            #Qt origin is at the top-left corner. The x values increase to the right and the y values increase downwards \
            #as compared to Genesis codinates where origin is center and y value is upwards, that is why ypos is negated 
            ypos = -(y-self.ymin)*self.yratio
        xpos = (x-self.xmin)*self.xratio
        
        return(xpos,ypos)

    def setupSlot(self,mooseObj,qgraphicItem):
        self.mooseId_GObj[element(mooseObj).getId()] = qgraphicItem
        qgraphicItem.connect(qgraphicItem,QtCore.SIGNAL("qgtextPositionChange(PyQt_PyObject)"),self.positionChange)
        qgraphicItem.connect(qgraphicItem,QtCore.SIGNAL("qgtextItemSelectedChange(PyQt_PyObject)"),self.objectEditSlot)

    def updateItemSlot(self, mooseObject):
        #This is overridden by derived classes to connect appropriate
        #slot for updating the display item.
        #In this case if the name is updated from the keyboard both in mooseobj and gui gets updation
        changedItem = ''

        for item in self.sceneContainer.items():
            #print "updateItemSlot",item
            if isinstance(item,PoolItem):
                if mooseObject.getId() == element(item.mobj).getId():
                    print "##",mooseObject.getId(),element(item.mobj).getId()
                    item.updateSlot()
                    #once the text is edited in editor, laydisplay gets updated in turn resize the length, positionChanged signal shd be emitted
                    self.positionChange(mooseObject)        
    def getToolBars(self):
        raise NotImplementedError('method must be reimplemented in subclass')
        '''if not hasattr(self, '_insertToolBar'):
            self._insertToolBar = QtGui.QToolBar('Insert')
            self._toolBars.append(self._insertToolBar)
            for action in self.insertMenu.actions():
                button = MToolButton()
                button.setDefaultAction(action)
                self._insertToolBar.addWidget(button)
        return self._toolBars
        '''
    def sizeHint(self):
        return QtCore.QSize(800,400)

    def updateModelView(self):
        # if self.modelRoot == '/':
        #     m = wildcardFind('/##[ISA=ChemCompt]')
        # else:
        #     m = wildcardFind(self.modelRoot+'/##[ISA=ChemCompt]')
        if not self.m:
            # when we want an empty GraphicView while creating new model,
            # then remove all the view and add an empty view
            if hasattr(self, 'view') and isinstance(self.view, QtGui.QWidget):
                self.layout().removeWidget(self.view)
            createdItem = {}
            self.sceneContainer.setSceneRect(-self.width()/2,-self.height()/2,self.width(),self.height())
            self.view = GraphicalView(self, self.modelRoot,self.sceneContainer,self.border,self,createdItem)
            self.connect(self.view, QtCore.SIGNAL("dropped"), self.objectEditSlot)
            hLayout = QtGui.QGridLayout(self)
            self.setLayout(hLayout)
            hLayout.addWidget(self.view)
        else:
            # maxmium and minimum coordinates of the objects specified in kkit file. 
            self.sceneContainer.clear()
            
            
            # Compartment and its members are put on the qgraphicsscene
            self.mooseObjOntoscene()
            
            # All the moose Object are connected for visualization 
            self.drawLine_arrow(itemignoreZooming=False)

            createdItem = {}
            if hasattr(self, 'view') and isinstance(self.view, QtGui.QWidget):
                self.layout().removeWidget(self.view)
            self.view = GraphicalView(self,self.modelRoot,self.sceneContainer,self.border,self,createdItem)
            hLayout = QtGui.QGridLayout(self)
            self.setLayout(hLayout)
            hLayout.addWidget(self.view)
            #self.layout().addWidget(self.view)
    
    '''
    def resetColor(self):
        for item in self.sceneContainer.items():
            if isinstance(item,PoolItem):
                pinfo = moose.element(item.mobj).path+'/info'
                color,bg = getColor(pinfo,self.colorMap)
                item.updateColor(bg)
    def changeBgSize(self):
        for item in self.sceneContainer.items():
            if isinstance(item,PoolItem):
                initialConc = moose.element(item.mobj).concInit
                presentConc = moose.element(item.mobj).conc
                if initialConc != 0:
                    ratio = presentConc/initialConc
                else:
                    # multipying by 1000 b'cos moose concentration is in milli in moose
                    ratio = presentConc*1000
                item.updateRect(math.sqrt(ratio))
        
    def colorChange(self):
        #While simulation is running pool color are increased or decreased as per concentration level 
        for item in self.sceneContainer.items():
            if isinstance(item,PoolItem):
                bg = item.returnColor()
                initialConc = moose.element(item.mobj[0]).concInit
                presentConc = moose.element(item.mobj[0]).conc
                if initialConc != 0:
                    ratio = presentConc/initialConc
                else:
                    # multipying by 1000 b'cos moose concentration is in milli in moose
                    ratio = presentConc*1000
                #alpha between0-255
                alpha = 128*ratio
                # Limiting alpha values
                if int(math.floor(alpha)) > 255:
                    alpha = 255
                elif int(math.floor(alpha))< 40:
                    alpha = 40
                #only alpha value is changed
                bg =QtGui.QColor(bg.red(),bg.green(),bg.blue(),alpha)
                item.updateColor(bg)
     '''   
    def positionChange(self,mooseObject):
        #If the item position changes, the corresponding arrow's are calculated
        if isinstance(element(mooseObject),CubeMesh):
            for k, v in self.qGraCompt.items():
                #mesh = mooseObject.path+'/mesh[0]'
                #print " mesh ",mesh,k.path
                #if k.path == mesh:
                for rectChilditem in v.childItems():
                    self.updateArrow(rectChilditem)
        else:
            mobj = self.mooseId_GObj[mooseObject.getId()]
            self.updateArrow(mobj)
            for k, v in self.qGraCompt.items():
                rectcompt = v.childrenBoundingRect()
                v.setRect(rectcompt.x()-10,rectcompt.y()-10,(rectcompt.width()+20),(rectcompt.height()+20))

    def emitItemtoEditor(self,mooseObject):
        #self.emit(QtCore.SIGNAL("itemPressed(PyQt_PyObject)"),mooseObject)
        self.editObject.emit(mooseObject.path)

    def drawLine_arrow(self, itemignoreZooming=False):
        for inn,out in self.srcdesConnection.items():
            # self.srcdesConnection is dictionary which contains key,value \
            #    key is Enzyme or Reaction  and value [[list of substrate],[list of product]] (tuple)
            #    key is FuncBase and value is [list of pool] (list)

            #src = self.mooseId_GObj[inn]
            if isinstance(out,tuple):
                if len(out[0])== 0:
                    print inn.className + ':' +inn.name+ " doesn't output message"
                else:
                    src = self.mooseId_GObj[inn.getId()]
                    for items in (items for items in out[0] ):
                        des = self.mooseId_GObj[element(items[0]).getId()]
                        self.lineCord(src,des,items,itemignoreZooming)
                if len(out[1]) == 0:
                    print inn.className + ':' +inn.name+ " doesn't output message"
                else:
                    for items in (items for items in out[1] ):
                        des = self.mooseId_GObj[element(items[0]).getId()]
                        self.lineCord(src,des,items,itemignoreZooming)
            elif isinstance(out,list):
                if len(out) == 0:
                    print "Func pool doesn't have sumtotal"
                else:
                    src = self.mooseId_GObj[element(inn).getId()]
                    for items in (items for items in out ):
                        des = self.mooseId_GObj[element(items[0]).getId()]
                        self.lineCord(src,des,items,itemignoreZooming)
    
    def lineCord(self,src,des,type_no,itemignoreZooming):
        endtype = type_no[1]
        line = 0
        if (src == "") and (des == ""):
            print "Source or destination is missing or incorrect"
            return 
        srcdes_list = [src,des,endtype,line]
        arrow = calcArrow(srcdes_list,itemignoreZooming,self.iconScale)
        self.drawLine(srcdes_list,arrow)

        while(type_no[2] > 1 and line <= (type_no[2]-1)):
            srcdes_list =[src,des,endtype,line]
            arrow = calcArrow(srcdes_list,itemignoreZooming,self.iconScale)
            self.drawLine(srcdes_list,arrow)
            line = line +1

        if type_no[2] > 5:
            print "Higher order reaction will not be displayed"

    def drawLine(self,srcdes_list,arrow):
        src = srcdes_list[0]
        des = srcdes_list[1]
        endtype = srcdes_list[2]
        line = srcdes_list[3]
        source = element(next((k for k,v in self.mooseId_GObj.items() if v == src), None))
        for l,v,o in self.object2line[src]:
            if v == des and o ==line:
                l.setPolygon(arrow)
                arrowPen = l.pen()
                arrowPenWidth = self.arrowsize*self.iconScale
                arrowPen.setColor(l.pen().color())
                arrowPen.setWidth(arrowPenWidth)
                l.setPen(arrowPen)
                return
        qgLineitem = self.sceneContainer.addPolygon(arrow)
        pen = QtGui.QPen(QtCore.Qt.green, 0, Qt.Qt.SolidLine, Qt.Qt.RoundCap, Qt.Qt.RoundJoin)
        pen.setWidth(self.arrowsize)
        #pen.setCosmetic(True)
        # Green is default color moose.ReacBase and derivatives - already set above
        if  isinstance(source, EnzBase):
            if ( (endtype == 's') or (endtype == 'p')):
                pen.setColor(QtCore.Qt.red)
            elif(endtype != 'cplx'):
                p1 = (next((k for k,v in self.mooseId_GObj.items() if v == src), None))
                pinfo = p1.path+'/info'
                color,bgcolor = getColor(pinfo,self.colorMap)
                pen.setColor(color)
        elif isinstance(source, moose.PoolBase):
            pen.setColor(QtCore.Qt.blue)
        elif isinstance(source,moose.StimulusTable):
            pen.setColor(QtCore.Qt.yellow)
        self.lineItem_dict[qgLineitem] = srcdes_list
        self.object2line[ src ].append( ( qgLineitem, des,line) )
        self.object2line[ des ].append( ( qgLineitem, src,line ) )
        qgLineitem.setPen(pen)

    def updateArrow(self,qGTextitem):
        #if there is no arrow to update then return
        if qGTextitem not in self.object2line:
            return
        listItem = self.object2line[qGTextitem]
        for ql, va,order in self.object2line[qGTextitem]:
            srcdes = []
            srcdes = self.lineItem_dict[ql]
            # Checking if src (srcdes[0]) or des (srcdes[1]) is ZombieEnz,
            # if yes then need to check if cplx is connected to any mooseObject, 
            # so that when Enzyme is moved, cplx connected arrow to other mooseObject(poolItem) should also be updated
            if( type(srcdes[0]) == EnzItem):
                self.cplxUpdatearrow(srcdes[0])
            elif( type(srcdes[1]) == EnzItem):
                self.cplxUpdatearrow(srcdes[1])
            
            # For calcArrow(src,des,endtype,itemignoreZooming) is to be provided
            arrow = calcArrow(srcdes,self.itemignoreZooming,self.iconScale)
            ql.setPolygon(arrow)
    
    def cplxUpdatearrow(self,srcdes):
        # srcdes which is 'EnzItem' from this,get ChildItems are retrived (b'cos cplx is child of zombieEnz)
        #And cplxItem is passed for updatearrow

        #Note: Here at this point enzItem has just one child which is cplxItem and childItems returns, PyQt4.QtGui.QGraphicsEllipseItem,CplxItem
        #Assuming CplxItem is always[1], but still check if not[0], if something changes in structure one need to keep an eye.
        if (srcdes.childItems()[1],CplxItem):
            self.updateArrow(srcdes.childItems()[1])
        else:
            self.updateArrow(srcdes.childItems()[0])
    
    def keyPressEvent(self,event):
        # key1 = event.key() # key event does not distinguish between capital and non-capital letters
        key = event.text().toAscii().toHex()
        if key ==  '41': # 'A' fits the view to iconScale factor
            itemignoreZooming = False
            self.updateItemTransformationMode(itemignoreZooming)
            self.view.fitInView(self.sceneContainer.itemsBoundingRect().x()-10,self.sceneContainer.itemsBoundingRect().y()-10,self.sceneContainer.itemsBoundingRect().width()+20,self.sceneContainer.itemsBoundingRect().height()+20,Qt.Qt.IgnoreAspectRatio)
            self.drawLine_arrow(itemignoreZooming=False)
            
        elif (key == '2e'): # '.' key, lower case for '>' zooms in 
            self.view.scale(1.1,1.1)

        elif (key == '2c'): # ',' key, lower case for '<' zooms in
            self.view.scale(1/1.1,1/1.1)

        elif (key == '3c'): # '<' key. zooms-in to iconScale factor
            self.iconScale *= 0.8
            self.updateScale( self.iconScale )

        elif (key == '3e'): # '>' key. zooms-out to iconScale factor
            self.iconScale *= 1.25
            self.updateScale( self.iconScale )
            
        elif (key == '61'):  # 'a' fits the view to initial value where iconscale=1
            self.iconScale = 1
            self.updateScale( self.iconScale )
            self.view.fitInView(self.sceneContainer.itemsBoundingRect().x()-10,self.sceneContainer.itemsBoundingRect().y()-10,self.sceneContainer.itemsBoundingRect().width()+20,self.sceneContainer.itemsBoundingRect().height()+20,Qt.Qt.IgnoreAspectRatio)
                   
    def updateItemTransformationMode(self, on):
        for v in self.sceneContainer.items():
            if( not isinstance(v,ComptItem)):
                #if ( isinstance(v, PoolItem) or isinstance(v, ReacItem) or isinstance(v, EnzItem) or isinstance(v, CplxItem) ):
                if isinstance(v,KineticsDisplayItem):
                    v.setFlag(QtGui.QGraphicsItem.ItemIgnoresTransformations, on)

    def updateScale( self, scale ):
        for item in self.sceneContainer.items():
            if isinstance(item,KineticsDisplayItem):
                item.refresh(scale)
                #iteminfo = item.mobj.path+'/info'
                #xpos,ypos = self.positioninfo(iteminfo)
                xpos = item.scenePos().x()
                ypos = item.scenePos().y()

                if isinstance(item,ReacItem) or isinstance(item,EnzItem) or isinstance(item,MMEnzItem):
                     item.setGeometry(xpos,ypos, 
                                     item.gobj.boundingRect().width(), 
                                     item.gobj.boundingRect().height())
                elif isinstance(item,CplxItem):
                    item.setGeometry(item.gobj.boundingRect().width()/2,item.gobj.boundingRect().height(), 
                                     item.gobj.boundingRect().width(), 
                                     item.gobj.boundingRect().height())
                elif isinstance(item,PoolItem):
                     item.setGeometry(xpos, ypos,item.gobj.boundingRect().width()
                                     +PoolItem.fontMetrics.width('  '), 
                                     item.gobj.boundingRect().height())
                     item.bg.setRect(0, 0, item.gobj.boundingRect().width()+PoolItem.fontMetrics.width('  '), item.gobj.boundingRect().height())

        self.drawLine_arrow(itemignoreZooming=False)
        for k, v in self.qGraCompt.items():
            rectcompt = v.childrenBoundingRect()
            comptPen = v.pen()
            comptWidth =  self.defaultComptsize*self.iconScale
            comptPen.setWidth(comptWidth)
            v.setPen(comptPen)
            v.setRect(rectcompt.x()-comptWidth,rectcompt.y()-comptWidth,(rectcompt.width()+2*comptWidth),(rectcompt.height()+2*comptWidth))

class kineticEditorWidget(KineticsWidget):
    def __init__(self,*args):
        KineticsWidget.__init__(self,*args)
        self.insertMenu = QtGui.QMenu('&Insert')
        self._menus.append(self.insertMenu)
        self.insertMapper = QtCore.QSignalMapper(self)
        classlist = ['CubeMesh','CylMesh','Pool','FuncPool','SumFunc','Reac','Enz','MMenz','StimulusTable','Table']
        insertMapper, actions = self.getInsertActions(classlist)
        for action in actions:
            self.insertMenu.addAction(action)        

    def makePoolItem(self, poolObj, qGraCompt):
        return PoolItem(poolObj, qGraCompt)
        
    def getToolBars(self):
        if not hasattr(self, '_insertToolBar'):
            self._insertToolBar = QtGui.QToolBar('Insert')
            self._toolBars.append(self._insertToolBar)
            for action in self.insertMenu.actions():
                button = MToolButton()
                button.setDefaultAction(action)
                self._insertToolBar.addWidget(button)
        return self._toolBars

class kineticRunWidget(KineticsWidget):
    def __init__(self,*args):
        KineticsWidget.__init__(self,*args)

    def makePoolItem(self, poolObj, qGraCompt):
        return PoolItemCircle(poolObj, qGraCompt)

    def getToolBars(self):
        return self._toolBars

    def changeBgSize(self):
        #print "here in chageBgSize"
        for item in self.sceneContainer.items():
            if isinstance(item,PoolItemCircle):
                initialConc = moose.element(item.mobj).concInit
                presentConc = moose.element(item.mobj).conc
                if initialConc != 0:
                    ratio = presentConc/initialConc
                else:
                    # multipying by 1000 b'cos moose concentration is in milli in moose
                    ratio = presentConc
                
                item.updateRect(math.sqrt(ratio))
    def resetColor(self):
        for item in self.sceneContainer.items():
            if isinstance(item,PoolItemCircle):
                item.returnEllispeSize()
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    size = QtCore.QSize(1024 ,768)
    #modelPath = 'Kholodenko'
    modelPath = 'acc27'
    #modelPath = 'acc8'
    #modelPath = '3ARECB'
    #modelPath = '3AreacB'
    #modelPath = '5AreacB'
    itemignoreZooming = False
    try:
        filepath = '../../Demos/Genesis_files/'+modelPath+'.g'
        filepath = '/home/harsha/genesis_files/gfile/'+modelPath+'.g'
        print filepath
        f = open(filepath, "r")
        loadModel(filepath,'/'+modelPath)
        
        #moose.le('/'+modelPath+'/kinetics')
        dt = KineticsWidget()
        dt.modelRoot ='/'+modelPath
        ''' Loading moose signalling model in python '''
        #execfile('/home/harsha/BuildQ/Demos/Genesis_files/scriptKineticModel.py')
        #dt.modelRoot = '/model'
        
        dt.updateModelView()
        dt.show()
  
    except  IOError, what:
      (errno, strerror) = what
      print "Error number",errno,"(%s)" %strerror
      sys.exit(0)
    sys.exit(app.exec_())
