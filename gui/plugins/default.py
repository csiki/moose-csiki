# default.py --- 
# 
# Filename: default.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Tue Nov 13 15:58:31 2012 (+0530)
# Version: 
# Last-Updated: Thu Jul 18 10:35:00 2013 (+0530)
#           By: subha
#     Update #: 2244
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# The default placeholder plugin for MOOSE
# 
# 

# Change log:
# 
# 
# 
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
# 
# 

# Code:

import sys
import config
import pickle
import os
from collections import defaultdict
import numpy as np
import re
from PyQt4 import QtGui, QtCore
from PyQt4.Qt import Qt
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import moose
from moose import utils
import mtree
from mtoolbutton import MToolButton
from msearch import SearchWidget
from checkcombobox import CheckComboBox

from mplugin import MoosePluginBase, EditorBase, EditorWidgetBase, PlotBase, RunBase
#from defaultToolPanel import DefaultToolPanel

class MoosePlugin(MoosePluginBase):
    """Default plugin for MOOSE GUI"""
    def __init__(self, root, mainwindow):
        MoosePluginBase.__init__(self, root, mainwindow)

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
            self.editorView = MooseEditorView(self)
            # self.editorView.getCentralWidget().editObject.connect(self.mainWindow.objectEditSlot)
            self.currentView = self.editorView
        return self.editorView

    def getPlotView(self):
        if not hasattr(self, 'plotView'):
            self.plotView = PlotView(self)
        return self.plotView

    def getRunView(self):
        if not hasattr(self, 'runView') or self.runView is None:
            self.runView = RunView(self)
        return self.runView

    def getMenus(self):
        """Create a custom set of menus."""
        return self._menus


class MooseEditorView(EditorBase):
    """Default editor.

    """
    def __init__(self, plugin):
        EditorBase.__init__(self, plugin)
        self.__initMenus()
        self.__initToolBars()

    def __initMenus(self):
        editMenu = QtGui.QMenu('&Edit')
        for menu in self.getCentralWidget().getMenus():
            editMenu.addMenu(menu)
        self._menus.append(editMenu)

    def __initToolBars(self):
        for toolbar in self.getCentralWidget().getToolBars():
            self._toolBars.append(toolbar)

    def getToolPanes(self):
        return super(MooseEditorView, self).getToolPanes()

    def getLibraryPane(self):
        return super(MooseEditorView, self).getLibraryPane()

    def getOperationsWidget(self):
        return super(MooseEditorView, self).getOperationsPane()

    def getCentralWidget(self):
        """Retrieve or initialize the central widget.

        Note that we call the widget's setModelRoot() function
        explicitly with the plugin's modelRoot as the argument. This
        enforces an update of the widget display with the current
        modelRoot.

        This function should be overridden by any derived class as it
        has the editor widget class hard coded into it.

        """
        if self._centralWidget is None:
            self._centralWidget = DefaultEditorWidget()
            if hasattr(self._centralWidget, 'init'):
                self._centralWidget.init()
            self._centralWidget.setModelRoot(self.plugin.modelRoot)
        return self._centralWidget


class MooseTreeEditor(mtree.MooseTreeWidget):
    """Subclass of MooseTreeWidget to implement drag and drop events. It
    creates an element under the drop location using the dropped mime
    data as class name.

    """
    def __init__(self, *args):
        mtree.MooseTreeWidget.__init__(self, *args)

    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat('text/plain'):
            event.acceptProposedAction()

    def dragMoveEvent(self, event):
        if event.mimeData().hasFormat('text/plain'):
            event.acceptProposedAction()

    def dropEvent(self, event):
        """Insert an element of the specified class in drop location"""
        if not event.mimeData().hasFormat('text/plain'):
            return
        pos = event.pos()
        item = self.itemAt(pos)
        print "dropEvent",item
        print item.mobj.path
        try:
            self.insertChildElement(item, str(event.mimeData().text()))
            event.acceptProposedAction()
        except NameError:
            return            


class DefaultEditorWidget(EditorWidgetBase):
    """Editor widget for default plugin. 
    
    Plugin-writers should code there own editor widgets derived from
    EditorWidgetBase.

    It adds a toolbar for inserting moose objects into the element
    tree. The toolbar contains MToolButtons for moose classes.

    Signals: editObject - inherited from EditorWidgetBase , emitted
    with currently selected element's path as argument. Should be
    connected to whatever slot is responsible for firing the object
    editor in top level.

    """
    def __init__(self, *args):
        EditorWidgetBase.__init__(self, *args)
        layout = QtGui.QHBoxLayout()
        self.setLayout(layout)
        self.tree = MooseTreeEditor()
        self.tree.setAcceptDrops(True)
        self.getTreeMenu()
        self.layout().addWidget(self.tree)        

    def getTreeMenu(self):
        try:
            return self.treeMenu
        except AttributeError:
            self.treeMenu = QtGui.QMenu()
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self.treeMenu.exec_)
        # Inserting a child element
        self.insertMenu = QtGui.QMenu('Insert')
        self._menus.append(self.insertMenu)
        self.treeMenu.addMenu(self.insertMenu)
        self.insertMapper = QtCore.QSignalMapper(self)
        ignored_bases = ['ZPool', 'Msg', 'Panel', 'SolverBase', 'none']
        ignored_classes = ['ZPool','ZReac','ZMMenz','ZEnz','CplxEnzBase']
        classlist = [ch[0].name for ch in moose.element('/classes').children
                     if (ch[0].baseClass not in ignored_bases)
                     and (ch[0].name not in (ignored_bases + ignored_classes))
                     and not ch[0].name.startswith('Zombie')
                     and not ch[0].name.endswith('Base')                     
                 ]
        insertMapper, actions = self.getInsertActions(classlist)
        for action in actions:
            self.insertMenu.addAction(action)
        self.connect(insertMapper, QtCore.SIGNAL('mapped(const QString&)'), self.tree.insertElementSlot)
        self.editAction = QtGui.QAction('Edit', self.treeMenu)
        self.editAction.triggered.connect(self.editCurrentObjectSlot)
        self.tree.elementInserted.connect(self.elementInsertedSlot)
        self.treeMenu.addAction(self.editAction)
        return self.treeMenu

    def updateModelView(self):
        self.tree.recreateTree(root=self.modelRoot)
        # if current in self.tree.odict:
        #     self.tree.setCurrentItem(current)

    def updateItemSlot(self, mobj):
        """This should be overridden by derived classes to connect appropriate
        slot for updating the display item.

        """
        self.tree.updateItemSlot(mobj)

    def editCurrentObjectSlot(self):
        """Emits an `editObject(str)` signal with moose element path of
        currently selected tree item as argument

        """
        mobj = self.tree.currentItem().mobj
        self.editObject.emit(mobj.path)

    def sizeHint(self):
        return QtCore.QSize(400, 300)

    def getToolBars(self):
        if not hasattr(self, '_insertToolBar'):
            self._insertToolBar = QtGui.QToolBar('Insert')
            for action in self.insertMenu.actions():
                button = MToolButton()
                button.setDefaultAction(action)
                self._insertToolBar.addWidget(button)
            self._toolBars.append(self._insertToolBar)
        return self._toolBars



############################################################
#
# View for running a simulation and runtime visualization
#
############################################################


from mplot import CanvasWidget

class RunView(RunBase):
    """A default runtime view implementation. This should be
    sufficient for most common usage.
    
    canvas: widget for plotting

    dataRoot: location of data tables

    """
    def __init__(self, *args, **kwargs):
        RunBase.__init__(self, *args, **kwargs)
        self.getCentralWidget()
        self.modelRoot = '/model'
        self.dataRoot = '/data'
        self.getCentralWidget()
        self.setModelRoot(moose.Neutral(self.plugin.modelRoot).path)
        self.setDataRoot(moose.Neutral('/data').path)
        self.plugin.modelRootChanged.connect(self.setModelRoot)
        self.plugin.dataRootChanged.connect(self.setDataRoot)        
        self._menus += self.getCentralWidget().getMenus()

    def getCentralWidget(self):
        """TODO: replace this with an option for multiple canvas
        tabs"""
        if self._centralWidget is None:
            self._centralWidget = PlotWidget()
        return self._centralWidget
        
    # def setDataRootSlot(self):
    #     path, ok = QtGui.QInputDialog.getText(self.getCentralWidget(), 'Set data root', 'Enter path to data root')
    #     if ok:
    #         self.setDataRoot(str(path))
        
    def setDataRoot(self, path):        
        self.dataRoot = path
        self.getCentralWidget().setDataRoot(path)
        self.getSchedulingDockWidget().widget().setDataRoot(path)

    def setModelRoot(self, path):
        self.modelRoot = path
        self.getSchedulingDockWidget().widget().setModelRoot(path)
        self.getCentralWidget().setModelRoot(path)

    def getDataTablesPane(self):
        """This should create a tree widget with dataRoot as the root
        to allow visual selection of data tables for plotting."""
        raise NotImplementedError()

    def plotAllData(self):
        """This is wrapper over the same function in PlotWidget."""
        self.centralWidget.plotAllData()

    def getToolPanes(self):
        if not self._toolPanes:
            self._toolPanes = [self.getSchedulingDockWidget()]
        return self._toolPanes

    def getSchedulingDockWidget(self):
        """Create and/or return a widget for schduling"""
        if hasattr(self, 'schedulingDockWidget')  and self.schedulingDockWidget is not None:
            return self.schedulingDockWidget
        self.schedulingDockWidget = QtGui.QDockWidget('Scheduling')
        widget = SchedulingWidget()
        widget.setDataRoot(self.dataRoot)
        widget.setModelRoot(self.modelRoot)
        self.schedulingDockWidget.setWidget(widget)
        widget.runner.update.connect(self._centralWidget.updatePlots)
        widget.runner.finished.connect(self._centralWidget.rescalePlots)
        widget.simtimeExtended.connect(self._centralWidget.extendXAxes)
        widget.runner.resetAndRun.connect(self._centralWidget.plotAllData)
        return self.schedulingDockWidget
    

class MooseRunner(QtCore.QObject):
    """Helper class to control simulation execution

    See: http://doc.qt.digia.com/qq/qq27-responsive-guis.html :
    'Solving a Problem Step by Step' for design details.
    """
    resetAndRun = QtCore.pyqtSignal(name='resetAndRun')
    update = QtCore.pyqtSignal(name='update')
    finished = QtCore.pyqtSignal(name='finished')
    '''
    _instance = None # singleton
    inited = False
    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(MooseRunner, cls).__new__(cls, *args, **kwargs)
        return cls._instance        
    '''    
    def __init__(self, *args, **kwargs):
        QtCore.QObject.__init__(self, *args, **kwargs)
        '''
        if (MooseRunner.inited):
            return
        '''
        self._updateInterval = 100e-3
        self._simtime = 0.0        
        self._clock = moose.Clock('/clock')
        self._pause = False
        self.dataRoot = moose.Neutral('/data').path
        self.modelRoot = moose.Neutral('/model').path
        #MooseRunner.inited = True

    def doResetAndRun(self, tickDtMap, tickTargetMap, simtime, updateInterval):
        self._pause = False
        self._updateInterval = updateInterval
        self._simtime = simtime
        utils.updateTicks(tickDtMap)
        utils.assignTicks(tickTargetMap)
        self.resetAndRun.emit()
        moose.reinit()
        QtCore.QTimer.singleShot(0, self.run)
        
    def run(self):
        """Run simulation for a small interval."""
        if self._clock.currentTime >= self._simtime:
            self.finished.emit()
            return
        if self._pause:
            return
        toRun = self._simtime - self._clock.currentTime
        if toRun > self._updateInterval:
            toRun = self._updateInterval
        moose.start(toRun)
        self.update.emit()
        QtCore.QTimer.singleShot(0, self.run)
    
    def continueRun(self, simtime, updateInterval):
        """Continue running without reset for `simtime`."""
        self._simtime = simtime
        self._updateInterval = updateInterval
        self._pause = False
        QtCore.QTimer.singleShot(0, self.run)
    
    def stop(self):
        """Pause simulation"""
        self._pause = True

class SchedulingWidget(QtGui.QWidget):
    """Widget for scheduling.

    Important member fields:

    runner - object to run/pause/continue simulation. Whenevr an
    `updateInterval` time has been simulated this object sends an
    `update()` signal. This can be connected to other objects to
    update their data.

    SIGNALS:
    resetAndRun(tickDt, tickTargets, simtime, updateInterval)

        tickDt: dict mapping tick nos to dt
        tickTargets: dict mapping ticks to target paths
        simtime: total simulation runtime
        updateInterval: interval between update signals are to be emitted.

    simtimeExtended(simtime)
        emitted when simulation time is increased by user.
    

    """    
    resetAndRun = QtCore.pyqtSignal(dict, dict, float, float, name='resetAndRun')
    simtimeExtended = QtCore.pyqtSignal(float, name='simtimeExtended')
    continueRun = QtCore.pyqtSignal(float, float, name='continueRun')
    def __init__(self, *args, **kwargs):
        QtGui.QWidget.__init__(self, *args, **kwargs)
        layout = QtGui.QVBoxLayout()
        self.advanceOptiondisplayed = False
        self.simtimeWidget = self.__getSimtimeWidget()
        self.tickListWidget = self.__getTickListWidget()
        self.runControlWidget = self.__getRunControlWidget()
        self.advanceOpt = self.__getAdvanceOptionsWidget()
        layout.addWidget(self.runControlWidget)
        layout.addWidget(self.simtimeWidget)
        layout.addWidget(self.advanceOpt)
        layout.addWidget(self.tickListWidget)

        if not self.advanceOptiondisplayed:
            self.tickListWidget.hide()

        self.updateInterval = 100e-3 # This will be made configurable with a textbox
        self.__getUpdateIntervalWidget()
        #layout.addWidget(self.__getUpdateIntervalWidget())
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        layout.addItem(spacerItem)
        self.setLayout(layout)
        self.runner = MooseRunner()
        self.resetAndRunButton.clicked.connect(self.resetAndRunSlot)
        self.continueButton.clicked.connect(self.doContinueRun)
        self.continueRun.connect(self.runner.continueRun)
        self.stopButton.clicked.connect(self.runner.stop)
        self.runner.update.connect(self.updateCurrentTime)
    def updateTickswidget(self):

        if self.advanceOptiondisplayed:
            self.tickListWidget.hide()
            self.advanceOptiondisplayed = False
        else:
            self.tickListWidget.show()
            self.advanceOptiondisplayed = True

    def __getAdvanceOptionsWidget(self):
        widget = QtGui.QWidget()
        layout = QtGui.QVBoxLayout()
        icon = QtGui.QIcon(os.path.join(config.settings[config.KEY_ICON_DIR],'arrow.png'))
        self.advancedOption = QtGui.QToolButton()
        self.advancedOption.setText("Advance Options")
        self.advancedOption.setIcon(QtGui.QIcon(icon))
        self.advancedOption.setToolButtonStyle( Qt.ToolButtonTextBesideIcon );
        self.advancedOption.clicked.connect(self.updateTickswidget)
        layout.addWidget(self.advancedOption)
        widget.setLayout(layout)
        return widget
    def __getUpdateIntervalWidget(self):
        label = QtGui.QLabel('Plot update interval')
        self.updateIntervalText = QtGui.QLineEdit(str(self.updateInterval))
        label = QtGui.QLabel('Update plots after every')
        ulabel = QtGui.QLabel('seconds of simulation')
        layout = QtGui.QHBoxLayout()
        layout.addWidget(label)
        layout.addWidget(self.updateIntervalText)
        layout.addWidget(ulabel)
        widget = QtGui.QWidget()
        widget.setLayout(layout)
        return widget
    
    def __getRunControlWidget(self):
        widget = QtGui.QWidget()
        layout = QtGui.QVBoxLayout()
        self.resetAndRunButton = QtGui.QPushButton('Reset and run')
        self.stopButton = QtGui.QPushButton('Stop')
        self.continueButton = QtGui.QPushButton('Continue')
        layout.addWidget(self.resetAndRunButton)
        layout.addWidget(self.stopButton)
        layout.addWidget(self.continueButton)
        widget.setLayout(layout)
        return widget

    def updateUpdateInterval(self):
        """Read the updateInterval from text box.

        If updateInterval is less than the smallest dt, then make it
        equal.
        """
        '''try:
            self.updateInterval = float(str(self.updateIntervalText.text()))
        except ValueError:
            QtGui.QMessageBox.warning(self, 'Invalid value', 'Specified plot update interval is meaningless.')
        '''
        #Harsha: Atleast for loading signalling model in the GSL method, the updateInterval need to be atleast
        #        equal to the min TickDt and not zero.
        tickDt = self.getTickDtMap().values()
        tickDt = [item for item in self.getTickDtMap().values() if float(item) != 0.0]
        dt = min(tickDt)
        #dt = min(self.getTickDtMap().values())
        if dt > self.updateInterval:
            self.updateInterval = dt
            #self.updateIntervalText.setText(str(dt))

    # def disableButton(self):
    #     """ When RunAndResetButton,continueButton,RunTillEndButton are clicked then disabling these buttons
    #     for further clicks"""
    #     self.disableButtons(False)

    # def enableButton(self):
    #     """ Enabling RunAndResetButton,continueButton,RunTillEndButton after stop button """
    #     self.disableButtons()

    # def disableButtons(self,Enabled=True):
    #     self.resetAndRunButton.setEnabled(Enabled)
    #     self.continueButton.setEnabled(Enabled)
    #     self.runTillEndButton.setEnabled(Enabled)

    def resetAndRunSlot(self):
        """This is just for adding the arguments for the function
        MooseRunner.resetAndRun"""
        self.updateUpdateInterval()
        tickDtMap = self.getTickDtMap()
        tickTargetsMap = self.getTickTargets()
        simtime = self.getSimTime()
        self.simtimeExtended.emit(simtime)
        self.runner.doResetAndRun(
            self.getTickDtMap(), 
            self.getTickTargets(), 
            self.getSimTime(), 
            self.updateInterval)

    def doContinueRun(self):
        """Helper function to emit signal with arguments"""
        
        #self.updateUpdateInterval()
        simtime = self.getSimTime()
        self.simtimeExtended.emit(simtime)
        self.continueRun.emit(simtime,
                               self.updateInterval)
    
    def __getSimtimeWidget(self):
        runtime = moose.Clock('/clock').runTime
        layout = QtGui.QGridLayout()
        simtimeWidget = QtGui.QWidget()
        self.simtimeEdit = QtGui.QLineEdit('1')
        self.simtimeEdit.setText(str(runtime))
        self.currentTimeLabel = QtGui.QLabel('0')
        layout.addWidget(QtGui.QLabel('Run for'), 0, 0)
        layout.addWidget(self.simtimeEdit, 0, 1)
        layout.addWidget(QtGui.QLabel('seconds'), 0, 2)        
        layout.addWidget(QtGui.QLabel('Current time:'), 1, 0)
        layout.addWidget(self.currentTimeLabel, 1, 1)
        layout.addWidget(QtGui.QLabel('second'), 1, 2)        
        simtimeWidget.setLayout(layout)
        return simtimeWidget

    def __getTickListWidget(self):
        layout = QtGui.QGridLayout()
        # Set up the column titles
        layout.addWidget(QtGui.QLabel('Tick'), 0, 0)
        layout.addWidget(QtGui.QLabel('dt'), 0, 1)
        layout.addWidget(QtGui.QLabel('Targets (wildcard)'), 0, 2, 1, 2)
        layout.setRowStretch(0, 1)
        # Create one row for each tick. Somehow ticks.shape is
        # (16,) while only 10 valid ticks exist. The following is a hack
        clock = moose.element('/clock')
        numticks = clock.numTicks
        for ii in range(numticks):
            tt = clock.tickDt[ii]
            layout.addWidget(QtGui.QLabel("(\'"+clock.path+'\').tickDt['+str(ii)+']'), ii+1, 0)
            layout.addWidget(QtGui.QLineEdit(str(tt)), ii+1, 1)
            layout.addWidget(QtGui.QLineEdit(''), ii+1, 2, 1, 2)
            layout.setRowStretch(ii+1, 1)
        '''
        for ii in range(ticks[0].localNumField):
            tt = ticks[ii]
            layout.addWidget(QtGui.QLabel(tt.path), ii+1, 0)
            layout.addWidget(QtGui.QLineEdit(str(tt.dt)), ii+1, 1)
            layout.addWidget(QtGui.QLineEdit(''), ii+1, 2, 1, 2)
            layout.setRowStretch(ii+1, 1)
        '''
        # We add spacer items to the last row so that expansion
        # happens at bottom. All other rows have stretch = 1, and
        # the last one has 0 (default) so that is the one that
        # grows
        rowcnt = layout.rowCount()
        for ii in range(3):
            layout.addItem(QtGui.QSpacerItem(1, 1), rowcnt, ii)
        layout.setRowStretch(rowcnt, 10)
        # layout.setColumnStretch(1, 1)
        layout.setColumnStretch(2, 2)
        widget = QtGui.QWidget()
        widget.setLayout(layout)
        return widget

    def updateCurrentTime(self):
        sys.stdout.flush()
        self.currentTimeLabel.setText('%f' % (moose.Clock('/clock').currentTime))

    def updateTextFromTick(self, tickNo):
        tick = moose.vector('/clock/tick')[tickNo]
        widget = self.tickListWidget.layout().itemAtPosition(tickNo + 1, 1).widget()
        if widget is not None and isinstance(widget, QtGui.QLineEdit):
            widget.setText(str(tick.dt))

    def updateFromMoose(self):
        """Update the tick dt from the tick objects"""
        ticks = moose.vector('/clock/tick')
        # Items at position 0 are the column headers, hence ii+1
        for ii in range(ticks[0].localNumField):
            self.updateTextFromTick(ii)
        self.updateCurrentTime()

    def getSimTime(self):
        try:
            time = float(str(self.simtimeEdit.text()))
            return time
        except ValueError, e:
            QtGui.QMessageBox.warning(self, 'Invalid value', 'Specified runtime was meaningless.')
        return 0

    def getTickTargets(self):
        """Return a dictionary containing tick nos mapped to the
        target specified in tick list widget. If target is empty, the
        tick is not included."""
        ret = {}
        for ii in range(1, self.tickListWidget.layout().rowCount()):
            widget = self.tickListWidget.layout().itemAtPosition(ii, 2).widget()
            if widget is not None and isinstance(widget, QtGui.QLineEdit):
                target = str(widget.text()).strip()
                if len(target) > 0:
                    ret[ii-1] = target
        print '77777', ret
        return ret

    def getTickDtMap(self):
        ret = {}
        # Items at position 0 are the column headers, hence ii+1
        for ii in range(1, self.tickListWidget.layout().rowCount()):
            widget = self.tickListWidget.layout().itemAtPosition(ii, 1).widget()
            if widget is not None and isinstance(widget, QtGui.QLineEdit):
                try:
                    ret[ii-1] = float(str(widget.text()))
                except ValueError:
                    QtGui.QMessageBox.warning(self, 'Invalid value', '`dt` for tick %d was meaningless.' % (ii-1))
        print '66666', ret
        return ret

    def setDataRoot(self, root='/data'):
        self.runner.dataRoot = moose.element(root).path

    def setModelRoot(self, root='/model'):
        self.runner.modelRoot = moose.element(root).path
        

from collections import namedtuple

# Keeps track of data sources for a plot. 'x' can be a table path or
# '/clock' to indicate time points from moose simulations (which will
# be created from currentTime field of the `/clock` element and the
# number of dat points in 'y'. 'y' should be a table. 'z' can be empty
# string or another table or something else. Will not be used most of
# the time (unless 3D or heatmap plotting).

PlotDataSource = namedtuple('PlotDataSource', ['x', 'y', 'z'], verbose=False)
event = None
legend = None
canvas = None
class PlotWidget(QtGui.QWidget):
    """A wrapper over CanvasWidget to handle additional MOOSE-specific
    stuff.

    modelRoot - path to the entire model our plugin is handling

    dataRoot - path to the container of data tables.

    TODO: do we really need this separation or should we go for
    standardizing location of data with respect to model root.

    pathToLine - map from moose path to Line2D objects in plot. Can
    one moose table be plotted multiple times? Maybe yes (e.g., when
    you want multiple other tables to be compared with the same data).

    lineToDataSource - map from Line2D objects to moose paths

    """
    def __init__(self, *args, **kwargs):
        QtGui.QWidget.__init__(self, *args)
        self.canvas = CanvasWidget()
        self.canvas.setParent(self)
        #global canvas
        #canvas = self.canvas
        self.navToolbar = NavigationToolbar(self.canvas, self)
        layout = QtGui.QGridLayout()
        layout.addWidget(self.canvas,0,1)
        layout.addWidget(self.navToolbar,1,1)
        self.setLayout(layout)
        self.modelRoot = '/'
        self.pathToLine = defaultdict(set)
        self.lineToDataSource = {}
        self.canvas.addSubplot(1, 1)
        self.onclick_count = 0
    @property
    def plotAll(self):
        return len(self.pathToLine) == 0

    def setModelRoot(self, path):
        self.modelRoot = path

    def setDataRoot(self, path):
        self.dataRoot = path

    def plotAllData(self):
        """Plot data from all tables under dataRoot"""        
        path = moose.element(self.dataRoot).path
        modelroot = moose.element(self.modelRoot).path
        time = moose.Clock('/clock').currentTime
        tabList = []
	#print " default ",path
        for tabId in moose.wildcardFind('%s/##[TYPE=Table]' % (path)):
            tab = moose.Table(tabId)
            line_list=[]
            tableObject = tab.neighbors['requestOut']
            if len(tableObject) > 0:
                # This is the default case: we do not plot the same
                # table twice. But in special cases we want to have
                # multiple variations of the same table on different
                # axes.
                #
                #Harsha: Adding color to graph for signalling model, check if given path has cubemesh or cylmesh
                color = ''
                if (len(moose.wildcardFind('%s/##[ISA=ChemCompt]' %(modelroot)))):
                    species = tableObject[0].path+'/info'
                    colormap_file = open(os.path.join(config.settings[config.KEY_COLORMAP_DIR], 'rainbow2.pkl'),'rb')
                    self.colorMap = pickle.load(colormap_file)
                    colormap_file.close()
                    hexchars = "0123456789ABCDEF"
                    #Genesis model exist the path and color will be set but not xml file so bypassing
                    if moose.exists(species):
                        color = moose.element(species).getField('color')
                        if ((not isinstance(color,(list,tuple)))):
                            if color.isdigit():
                                tc = int(color)
                                tc = (tc * 2 )
                                r,g,b = self.colorMap[tc]
                                color = "#"+ hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]
                lines = self.pathToLine[tab.path]

                if len(lines) == 0:
                    #Harsha: pass color for plot if exist and not white else random color
                    if (color != 'white'):
                        newLines = self.addTimeSeries(tab, label=tab.name,color=color)
                    else:
                        newLines = self.addTimeSeries(tab, label=tab.name)
                    self.pathToLine[tab.path].update(newLines)
                    for line in newLines:
                        self.lineToDataSource[line] = PlotDataSource(x='/clock', y=tab.path, z='')
                else:
                    for line in lines:
                        dataSrc = self.lineToDataSource[line]
                        xSrc = moose.element(dataSrc.x)
                        ySrc = moose.element(dataSrc.y)
                        if isinstance(xSrc, moose.Clock):
                            ts = np.linspace(0, time, len(tab.vector))
                        elif isinstance(xSrc, moose.Table):
                            ts = xSrc.vector.copy()
                        line.set_data(ts, tab.vector.copy())
                tabList.append(tab)
        self.canvas.mpl_connect('pick_event',self.onclick)
        if len(tabList) > 0:
            leg = self.canvas.callAxesFn('legend',loc='upper center',prop={'size':10},bbox_to_anchor=(0.5, -0.03),fancybox=True, shadow=True, ncol=3)
            #leg = self.canvas.callAxesFn('legend')
            #leg = self.canvas.callAxesFn('legend',loc='upper left', fancybox=True, shadow=True)
            global legend
            legend =leg
            for legobj in leg.legendHandles:
                legobj.set_linewidth(4.0)
                legobj.set_picker(True)
        self.canvas.draw()
    
    def onclick(self,event1):
        #print "onclick",event1.artist.get_label()
        #harsha:To workout with double-event-registered on onclick event
        #http://stackoverflow.com/questions/16278358/double-event-registered-on-mouse-click-if-legend-is-outside-axes
        if self.onclick_count % 2 == 0:
            legline = event1.artist
            #vis = event1.artist.get_visible()
            #self.canvas.figure.get_axes()[0].lines[4].set_visible(True)
            axes = self.canvas.figure.get_axes()
            for a in range(len(axes)):
                #lines =self.canvas.figure.get_axes()[a].lines
                lines = axes[a].lines
                for plotline in lines:
                    if plotline.get_label() == event1.artist.get_label():
                        vis = not plotline.get_visible()
                        plotline.set_visible(vis)
            #global event
            #event = event1
            if vis:
                legline.set_alpha(1.0)
            else:
                legline.set_alpha(0.2)
            self.canvas.draw()
        self.onclick_count+=1

        '''leg = self.canvas.callAxesFn('legend',loc='upper center',prop={'size':10},bbox_to_anchor=(0.5, -0.03),fancybox=True, shadow=True, ncol=3)
        print dir(leg)
        for l in leg.get_lines():
            l.set_visible(vis)
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        '''
    def addTimeSeries(self, table, *args, **kwargs):        
        ts = np.linspace(0, moose.Clock('/clock').currentTime, len(table.vector))
        return self.canvas.plot(ts, table.vector, *args, **kwargs)
        
    def addRasterPlot(self, eventtable, yoffset=0, *args, **kwargs):
        """Add raster plot of events in eventtable.

        yoffset - offset along Y-axis.
        """
        y = np.ones(len(eventtable.vector)) * yoffset
        return self.canvas.plot(eventtable.vector, y, '|')

    def updatePlots(self):
        for path, lines in self.pathToLine.items():            
            tab = moose.Table(path)
            data = tab.vector
            ts = np.linspace(0, moose.Clock('/clock').currentTime, len(data))
            for line in lines:
                line.set_data(ts, data)
        self.canvas.draw()

    def extendXAxes(self, xlim):
        for axes in self.canvas.axes.values():
            # axes.autoscale(False, axis='x', tight=True)
            axes.set_xlim(right=xlim)
            axes.autoscale_view(tight=True, scalex=True, scaley=True)
        self.canvas.draw()

    def rescalePlots(self):
        """This is to rescale plots at the end of simulation.
        
        ideally we should set xlim from simtime.
        """
        for axes in self.canvas.axes.values():
            axes.autoscale(True, tight=True)
            axes.relim()
            axes.autoscale_view(tight=True,scalex=True,scaley=True)
        self.canvas.draw()
    #Harsha: Passing directory path to save plots
    def saveCsv(self, line,directory):
        """Save selected plot data in CSV file"""
        src = self.lineToDataSource[line]
        xSrc = moose.element(src.x)
        ySrc = moose.element(src.y)
        y = ySrc.vector.copy()
        if isinstance(xSrc, moose.Clock):
            x = np.linspace(0, xSrc.currentTime, len(y))
        elif isinstance(xSrc, moose.Table):
            x = xSrc.vector.copy()
        filename = str(directory)+'/'+'%s.csv' % (ySrc.name)
        np.savetxt(filename, np.vstack((x, y)).transpose())
        print 'Saved data from %s and %s in %s' % (xSrc.path, ySrc.path, filename)

    def saveAllCsv(self):
        """Save data for all currently plotted lines"""
        #Harsha: Plots were saved in GUI folder instead provided QFileDialog box to save to
        #user choose
        fileDialog2 = QtGui.QFileDialog(self)
        fileDialog2.setFileMode(QtGui.QFileDialog.Directory)
        fileDialog2.setWindowTitle('Select Directory to save plots')
        fileDialog2.setOptions(QtGui.QFileDialog.ShowDirsOnly)
        fileDialog2.setLabelText(QtGui.QFileDialog.Accept, self.tr("Save"))
        targetPanel = QtGui.QFrame(fileDialog2)
        targetPanel.setLayout(QtGui.QVBoxLayout())
        layout = fileDialog2.layout()
        layout.addWidget(targetPanel)
        if fileDialog2.exec_():
            directory = fileDialog2.directory().path()
            for line in self.lineToDataSource.keys():
                self.saveCsv(line,directory)

    def getMenus(self):
        if not hasattr(self, '_menus'):
            self._menus = []
            self.plotAllAction = QtGui.QAction('Plot all data', self)
            self.plotAllAction.triggered.connect(self.plotAllData)
            self.plotMenu = QtGui.QMenu('Plot')
            self.plotMenu.addAction(self.plotAllAction)
            self.saveAllCsvAction = QtGui.QAction('Save all data in CSV files', self)
            self.saveAllCsvAction.triggered.connect(self.saveAllCsv)
            self.plotMenu.addAction(self.saveAllCsvAction)
            self._menus.append(self.plotMenu)
        return self._menus


###################################################
#
# Plot view - select fields to record
#
###################################################

class PlotView(PlotBase):
    """View for selecting fields on elements to plot."""
    def __init__(self, *args):
        PlotBase.__init__(self, *args)
        self.plugin.modelRootChanged.connect(self.getSelectionPane().setSearchRoot)
        self.plugin.dataRootChanged.connect(self.setDataRoot)
        self._recordingDict = {}
        self._reverseDict = {}
        self.dataRoot = self.plugin.dataRoot

    def setDataRoot(self, root):
        self.dataRoot = moose.element(root).path

    def getToolPanes(self):
        return (self.getFieldSelectionDock(), )

    def getSelectionPane(self):
        """Creates a widget to select elements and fields for plotting.
        search-root, field-name, comparison operator , value
        """
        if not hasattr(self, '_selectionPane'):
            self._searchWidget = SearchWidget()
            self._searchWidget.setSearchRoot(self.plugin.modelRoot)
            self._fieldLabel = QtGui.QLabel('Field to plot')
            self._fieldEdit = QtGui.QLineEdit()
            self._fieldEdit.returnPressed.connect(self._searchWidget.searchSlot)
            self._selectionPane = QtGui.QWidget()
            layout = QtGui.QHBoxLayout()
            layout.addWidget(self._fieldLabel)
            layout.addWidget(self._fieldEdit)
            self._searchWidget.layout().addLayout(layout)
            self._selectionPane = self._searchWidget
            self._selectionPane.layout().addStretch(1)
        return self._selectionPane

    def getOperationsPane(self):
        """TODO: complete this"""
        if hasattr(self, 'operationsPane'):
            return self.operationsPane
        self.operationsPane = QtGui.QWidget()
        self._createTablesButton = QtGui.QPushButton('Create tables for recording selected fields', self.operationsPane)
        self._createTablesButton.clicked.connect(self.setupRecording)
        layout = QtGui.QVBoxLayout()
        self.operationsPane.setLayout(layout)
        layout.addWidget(self._createTablesButton)
        return self.operationsPane

    def getFieldSelectionDock(self):
        if not hasattr(self, '_fieldSelectionDock'):
            self._fieldSelectionDock = QtGui.QDockWidget('Search and select elements')
            self._fieldSelectionWidget = QtGui.QWidget()
            layout = QtGui.QVBoxLayout()
            self._fieldSelectionWidget.setLayout(layout)
            layout.addWidget(self.getSelectionPane())
            layout.addWidget(self.getOperationsPane())
            self._fieldSelectionDock.setWidget(self._fieldSelectionWidget)
        return self._fieldSelectionDock

    def getCentralWidget(self):
        if not hasattr(self, '_centralWidget') or self._centralWidget is None:
            self._centralWidget = PlotSelectionWidget()
            self.getSelectionPane().executed.connect(self.selectElements)
        return self._centralWidget

    def selectElements(self, elements):
        """Refines the selection.
        
        Currently checks if _fieldEdit has an entry and if so, selects
        only elements which have that field, and ticks the same in the
        PlotSelectionWidget.

        """
        field = str(self._fieldEdit.text()).strip()
        if len(field) == 0:
            self.getCentralWidget().setSelectedElements(elements)
            return
        classElementDict = defaultdict(list)        
        for epath in elements:
            el = moose.element(epath)
            classElementDict[el.className].append(el)
        refinedList = []
        elementFieldList = []
        for className, elist in classElementDict.items():
            if field in elist[0].getFieldNames('valueFinfo'):
                refinedList +=elist
                elementFieldList += [(el, field) for el in elist]
        self.getCentralWidget().setSelectedElements(refinedList)
        self.getCentralWidget().setSelectedFields(elementFieldList)
        

    def setupRecording(self):
        """Create the tables for recording selected data and connect them."""
        for element, field in self.getCentralWidget().getSelectedFields():
            self.createRecordingTable(element, field)


    def createRecordingTable(self, element, field):
        """Create table to record `field` from element `element`

        Tables are created under `dataRoot`, the names are generally
        created by removing `/model` in the beginning of `elementPath`
        and replacing `/` with `_`. If this conflicts with an existing
        table, the id value of the target element (elementPath) is
        appended to the name.

        """
        if len(field) == 0 or ((element, field) in self._recordingDict):            
            return
        # The table path is not foolproof - conflict is
        # possible: e.g. /model/test_object and
        # /model/test/object will map to same table. So we
        # check for existing table without element field
        # path in recording dict.
        relativePath = element.path.partition('/model[0]/')[-1]
        if relativePath.startswith('/'):
            relativePath = relativePath[1:]
        #Convert to camelcase
        if field == "concInit":
            field = "ConcInit"
        elif field == "conc":
            field = "Conc"
        elif field == "nInit":
            field = "NInit"
        elif field == "n":
            field = "N"
        elif field == "volume":
            field = "Volume"
        elif field == "diffConst":
            field ="DiffConst"

        tablePath =  relativePath.replace('/', '_') + '.' + field
        tablePath = re.sub('.', lambda m: {'[':'_', ']':'_'}.get(m.group(), m.group()),tablePath)
        tablePath = self.dataRoot + '/' +tablePath
        if moose.exists(tablePath):
            tablePath = '%s_%d' % (tablePath, element.getId().value)
        if not moose.exists(tablePath):
            table = moose.Table(tablePath)
            print 'Created', table.path, 'for plotting', '%s.%s' % (element.path, field)
            target = element
            moose.connect(table, 'requestOut', target, 'get%s' % (field))
            self._recordingDict[(target, field)] = table
            self._reverseDict[table] = (target, field)

class PlotSelectionWidget(QtGui.QScrollArea):
    """Widget showing the fields of specified elements and their plottable
    fields. User can select any number of fields for plotting and click a
    button to generate the tables for recording data.

    The data tables are by default created under /data. One can call
    setDataRoot with a path to specify alternate location.

    """
    def __init__(self, *args):
        QtGui.QScrollArea.__init__(self, *args)        
        model = moose.Neutral('/model')
        self.modelRoot = model.path
        self.setLayout(QtGui.QVBoxLayout(self))
        self.layout().addWidget(self.getPlotListWidget())
        self.setDataRoot(moose.Neutral('/data'))
        self._elementWidgetsDict = {} # element path to corresponding qlabel and fields combo

    def getPlotListWidget(self):
        """An internal widget to display the list of elements and their
        plottable fields in comboboxes."""
        if not hasattr(self, '_plotListWidget'):
            self._plotListWidget = QtGui.QWidget(self)
            layout = QtGui.QGridLayout(self._plotListWidget)
            self._plotListWidget.setLayout(layout)        
            layout.addWidget(QtGui.QLabel('<h1>Elements matching search criterion will be listed here</h1>'), 0, 0)            
        return self._plotListWidget            
    
    def setSelectedElements(self, elementlist):
        """Create a grid of widgets displaying paths of elements in
        `elementlist` if it has at least one plottable field (a field
        with a numeric value). The numeric fields are listed in a
        combobox next to the element path and can be selected for
        plotting by the user.

        """
        for ii in range(self.getPlotListWidget().layout().count()):
            item = self.getPlotListWidget().layout().itemAt(ii)
            if item is None:
                continue
            self.getPlotListWidget().layout().removeItem(item)            
            w = item.widget()
            w.hide()
            del w
            del item
        self._elementWidgetsDict.clear()
        label = QtGui.QLabel('Element')
        label.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.getPlotListWidget().layout().addWidget(label, 0, 0, 1, 2)
        self.getPlotListWidget().layout().addWidget(QtGui.QLabel('Fields to plot'), 0, 2, 1, 1)
        for ii, entry in enumerate(elementlist):
            el = moose.element(entry)
            plottableFields = []
            for field, dtype in  moose.getFieldDict(el.className, 'valueFinfo').items():
                if dtype == 'double':
                    plottableFields.append(field)
            if len(plottableFields) == 0:
                continue
            elementLabel = QtGui.QLabel(el.path)
            fieldsCombo = CheckComboBox(self)
            fieldsCombo.addItem('')
            for item in plottableFields:
                fieldsCombo.addItem(item)
            self.getPlotListWidget().layout().addWidget(elementLabel, ii+1, 0, 1, 2)
            self.getPlotListWidget().layout().addWidget(fieldsCombo, ii+1, 2, 1, 1)
            self._elementWidgetsDict[el] = (elementLabel, fieldsCombo)
                
    def setModelRoot(self, root):
        pass

    def setDataRoot(self, path):
        """The tables will be created under dataRoot"""
        self.dataRoot = path

    def getSelectedFields(self):
        """Returns a list containing (element, field) for all selected fields"""
        ret = []
        for el, widgets in self._elementWidgetsDict.items():
            combo = widgets[1]
            for ii in range(combo.count()):
                field = str(combo.itemText(ii)).strip()
                if len(field) == 0:
                    continue
                checked, success = combo.itemData(ii, Qt.CheckStateRole).toInt()
                if success and checked == Qt.Checked:
                    ret.append((el, field))
        return ret

    def setSelectedFields(self, elementFieldList):
        """Set the checked fields for each element in elementFieldList.

        elementFieldList: ((element1, field1), (element2, field2), ...)
        
        """
        for el, field in elementFieldList:
            combo = self._elementWidgetsDict[el][1]
            idx = combo.findText(field)
            if idx >= 0:
                combo.setItemData(idx, QtCore.QVariant(Qt.Checked), Qt.CheckStateRole)
                combo.setCurrentIndex(idx)
# 
# default.py ends here
