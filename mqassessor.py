import tkinter as tk
from tkinter import ttk 
from tkinter import Toplevel
from tkinter import messagebox
from tkinter.filedialog import askopenfilename
import pandas as pd
import numpy as np
from math import log
from math import sqrt
from math import ceil
from scipy.optimize import minimize
from scipy.stats import chi2
from scipy.stats import t


def GUI():
    win = tk.Tk()
    win.title("MQAssessor")
    
    # creating a menu bar
    menu_bar = tk.Menu(win)
    win.config(menu=menu_bar)
    RWidth = 1024
    RHeight = 600
    win.geometry(("%dx%d")%(RWidth, RHeight))
    
    global isInputDone
    global isOutputDone
    global isModelDone
    isInputDone = False
    isOutputDone = False
    isModelDone = False

    file_menu = tk.Menu(menu_bar, tearoff=0)
    file_menu.add_command(label="Open raw data",
                          command=open_raw)
    file_menu.add_command(label="Open covariance data",
                          command=open_cov)
    file_menu.add_command(label="Open analysis",
                          command=open_model)
    
    #file_menu.add_separator()
    #file_menu.add_command(label="Save output as a txt file", command=saveoutput)
    file_menu.add_separator()
    file_menu.add_command(label="Save analysis",
                          command=isSaveOk)
    file_menu.add_command(label="Save output as a txt file",
                          command=isOutputOk)
    file_menu.add_command(label="Save output as an Excel file (only tables)",
                          command=isOutputOk2)
    file_menu.add_separator()
    file_menu.add_command(label="Exit", command=on_closing)
    menu_bar.add_cascade(label="File", menu=file_menu)
    
    var_menu = tk.Menu(menu_bar, tearoff=0)
    var_menu.add_command(label="Determine measurement model", command=isVarAssgn)
    menu_bar.add_cascade(label="Model", menu=var_menu)
    
    run_menu = tk.Menu(menu_bar, tearoff=0)
    run_menu.add_command(label="Run options", command=isRun)
    menu_bar.add_cascade(label="Run", menu=run_menu)
    
    help_menu = tk.Menu(menu_bar, tearoff=0)
    help_menu.add_command(label="About MQAssessor", command=about_mq)
    menu_bar.add_cascade(label="Help", menu=help_menu)
    
    #win.iconbitmap('taichi.ico')
    win.mainloop()

def open_raw():
    global isModelDone
    FileName = askopenfilename(initialdir="/", filetypes=(
            ["Excel File", ("*.xlsx", "*.xls")],
            ["csv file", "*.csv"]
            ))
    if "xls" in FileName:
        OpenFile = pd.read_excel(FileName)
    elif "csv" in FileName:
        OpenFile = pd.read_csv(FileName)
    Labels = OpenFile.columns
    Data = OpenFile.values
    obsvar(Labels, Data, isRaw=True)
    isModelDone = False


def open_cov():
    global isModelDone
    FileName = askopenfilename(initialdir="/", filetypes=(
            ["Excel File", ("*.xlsx", "*.xls")],
            ["csv file", "*.csv"]
            ))
    if "xls" in FileName:
        OpenFile = pd.read_excel(FileName)
    elif "csv" in FileName:
        OpenFile = pd.read_csv(FileName)
    Labels = OpenFile.columns
    Data = OpenFile.values
    obsvar(Labels, Data, isRaw=False)
    isModelDone = False


def open_model():
    FileName = askopenfilename(initialdir="/", filetypes=[(
        "MQAssessor File", "*.mqa")])
    OpenFile = pd.read_csv(FileName)
    Labels = list(OpenFile.columns)
    Data = list(OpenFile.values)
    global ObsNamesTemp
    global LtnNames
    global NumOfObsVarTemp
    global NumOfLtnVar
    global ObsCovsTemp
    global SampleSize
    global ObsToLtnTemp
    global isInputDone
    global isModelDone
    # Number of observed variables
    NumOfObsVarTemp = len(Labels) - 2
    ObsNamesTemp = NumOfObsVarTemp * [""]
    # Number of latent variables
    NumOfLtnVar = len(Data) - NumOfObsVarTemp - 1
    ObsCovsTemp = [[0] * NumOfObsVarTemp for i in range(NumOfObsVarTemp)]
    # Name of observed variables
    for i in range(NumOfObsVarTemp):
        ObsNamesTemp[i] = Labels[i + 2]
    # Observed covariances
        for j in range(NumOfObsVarTemp):
            ObsCovsTemp[i][j] = Data[i][j + 2]
    
    SampleSize = round(Data[NumOfObsVarTemp][2])
    ObsToLtnTemp = NumOfObsVarTemp * [0]
    LtnNames = NumOfLtnVar * [""]
    
    for i in range(NumOfLtnVar):
        LtnNames[i] = Data[NumOfObsVarTemp + 1 + i][1]
        for j in range(NumOfObsVarTemp):
            if round(Data[NumOfObsVarTemp + 1 + i][2+ j]) == 1:
                ObsToLtnTemp[j] = i
    isInputDone = True
    isModelDone = True
    VarAssgn()

def isSaveOk():
    if not isInputDone or NumOfLtnVar == 0:
        messagebox.showwarning("No data or latent variables", 
                               "Please input data or specify latent variables first.")
    else:
        save_model()

def isOutputOk():
    if not isOutputDone:
        messagebox.showwarning("No output", 
                               "Please get the output first.")
    else:
        save_output()

def isOutputOk2():
    if not isOutputDone:
        messagebox.showwarning("No output", 
                               "Please get the output first.")
    else:
        save_excel()
        
def save_model():
    # First, making a matrix
    #Mtr = np.empty([NumOfObsVarTemp + NumOfLtnVar + 2, NumOfObsVarTemp + 1])
    Mtr = [[""] * (NumOfObsVarTemp + 1) for i in range(NumOfObsVarTemp + NumOfLtnVar + 2)]
    Mtr[0][0] = " "
    for i in range(NumOfObsVarTemp):
        # The first row of the matrix is the name of the Observed variable
        Mtr[0][i + 1] = Mtr[i + 1][0] = ObsNames[i]
        # Next, The correlation matrix between the observed variables 
        for j in range(NumOfObsVarTemp):
            Mtr[i + 1][j + 1] = ObsCovsTemp[i][j]
    # Next, the sample size is added one line
    Mtr[NumOfObsVarTemp + 1][0] = "SampleSize"
    for i in range(NumOfObsVarTemp):
        Mtr[NumOfObsVarTemp + 1][i + 1] = SampleSize
    # Name of latent variables and its connection to observed variables
    for j in range(NumOfObsVarTemp):
        for i in range(NumOfLtnVar):
            Mtr[NumOfObsVarTemp + 2 + i][0] = LtnNames[i]
            if ObsToLtnTemp[j] == i:
                Mtr[NumOfObsVarTemp + 2 + i][j + 1] = 1
            else:
                Mtr[NumOfObsVarTemp + 2 + i][j + 1] = 0
    # converting to pandas dataframe
    Df = pd.DataFrame(Mtr)
    pd.DataFrame.to_csv(Df, header=False)
    filename = tk.filedialog.asksaveasfilename(defaultextension='.mqa',
                                               filetypes = [('MQAssessor files', '.mqa')])
    f = open(filename, 'w')
    a = pd.DataFrame.to_csv(Df, header=False)
    f.write(a)
    f.close()
        
def save_output():
    filename = tk.filedialog.asksaveasfilename(defaultextension='.txt',
                                                    filetypes = [('all files', '.*'), ('text files', '.txt')])
    f = open(filename, 'w')
    f.write(AllOutput)
    f.close()

def save_excel()    :
    filename = tk.filedialog.asksaveasfilename(defaultextension='.xlsx',
                                                    filetypes = [('Excel files', '.xlsx')])
    writer = pd.ExcelWriter(filename, engine = 'xlsxwriter')
    
    DfCoeffs = pd.DataFrame.from_dict(Coeffs)
    DfCoeffs.to_excel(writer, sheet_name="Lambda")
    DfResids = pd.DataFrame.from_dict(Resids)
    DfResids.to_excel(writer, sheet_name="Residual")
    DfCors = pd.DataFrame.from_dict(Cors)
    DfCors.to_excel(writer, sheet_name="Correlation")
    DfDVs = pd.DataFrame.from_dict(DVs)
    DfDVs.to_excel(writer, sheet_name="DiscriminantValidity")
    DfAVEs = pd.DataFrame.from_dict(AVEs)
    DfAVEs.to_excel(writer, sheet_name="AVE")
    DfDCs = pd.DataFrame.from_dict(DCs)
    DfDCs.to_excel(writer, sheet_name="DC")
    DfFPs = pd.DataFrame.from_dict(FPs)
    DfFPs.to_excel(writer, sheet_name="Pattern")
    DfFSs = pd.DataFrame.from_dict(FSs)
    DfFSs.to_excel(writer, sheet_name="Structure")
    DfRels = pd.DataFrame.from_dict(Rels)
    DfRels.to_excel(writer, sheet_name="Reliability")
    DfTbls = pd.DataFrame.from_dict(Tbls)
    DfTbls.to_excel(writer, sheet_name="CorrelationTable")
    writer.save()
    
    
def obsvar(Labels, Data, isRaw):
    # The ObsNamesTemp and NumOfObsVarTemp variables are used globally because they
    # are common to many functions.
    global ObsNamesTemp
    global NumOfObsVarTemp
    global ObsCovsTemp
    global SampleSize
    global isInputDone

    ObsNamesTemp = list(Labels)
    NumOfObsVarTemp = len(ObsNamesTemp)

    # A warning message is sent if the first line contains
    # only numeric elements.
    warning = False
    for i in range(NumOfObsVarTemp):
        if isNumber(ObsNamesTemp[i]):
            warning = True

    if warning:
        messagebox.showwarning("No variable name",
                               "The first line of the data should contain the name of each observed variable, not just numbers")

    # If the covariance matrix is unsymmetric about the diagonal line, issue
    # a warning message. If it is symmetric, it is confirmed as input data.
    if isRaw:
        # All elements in the source data must be numbers
        # I want to check it through isNumeric
        isNumeric = True
        # The data should not have a missing value. I want to check it out.
        isNoMissing = True
        for i in range(len(Data)):
            for j in range(len(Data[0])):
                if not isNumber(Data[i][j]):
                    isNumeric = False
                if pd.isnull(Data[i][j]):
                    isNoMissing = False

        if not isNumeric:
            messagebox.showwarning("Not numeric",
                               "Except for the first row, the data must be purely numeric.")
        elif not isNoMissing:
            messagebox.showwarning("Missing values",
                               "Missing values are not allowed. There should be a numerical value for all observed variables for all respondents.")
        else:
            SampleSize = len(Data)
            ObsCovsTemp = np.cov(np.transpose(Data))
            isInputDone = True
            
    # If the data is a covariance matrix
    else:
        isSymmetric = True
        for i in range(NumOfObsVarTemp):
            for j in range(i + 1, NumOfObsVarTemp):
                if Data[i][j] != Data[j][i]:
                    isSymmetric = False
        if isSymmetric:
            ObsCovsTemp = Data
            # Temporarily set it to -1 until samplesize is typed.
            SampleSize = -1
            isInputDone = True
            GetSampleSize()
        else:
            messagebox.showwarning("Not a covariance matrix",
                               "The covariance matrix should be symmetric above and below the diagonal, but this input is not.")

def about_mq():
    tk.messagebox.showinfo("About MQAssessor", "MQAsessor ver 0.0 \n" +
                           "Use of this program is free, subject to the following citation:\n" +
                           "Discriminant validity: What it is and how to assess it")


def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def on_closing():
    global win
    if messagebox.askokcancel("Quit", "Do you want to quit?"):
        win.destroy()


def isVarAssgn():
    if not isInputDone:
        messagebox.showwarning("No data", "Load the data file first.")
    else:
        VarAssgn()

def isRun():
    if not isModelDone:
        messagebox.showwarning("Model undetermined",
                               "Determine the measurement model first.")
    else:
        RunOption()

class GetSampleSize(tk.Tk):
    def __init__(self):

        tk.Tk.__init__(self)
        self.title("Sample size input")
        self.label = tk.Label(self,text = "Enter the sample size or the number of observations as an integer: ",
                              padx = 5, pady = 5)
        self.button = tk.Button(self, text = "OK", command = self.getme,
                                width = 5)
        
        self.entry = tk.Entry(self, width = 8) 
        self.label.grid(row=0, column = 0, padx = 10, pady = 10)
        self.entry.grid(row=0, column = 1, padx = 10, pady = 10)
        self.button.grid(row=0, column = 2, padx = 10, pady = 10)
        self.entry.bind("<Return>", self.getme)
        self.entry.focus()
        
    def getme(self, event=None):
        global SampleSize
        SampleSize = int(self.entry.get())
        self.destroy()


class VarAssgn(Toplevel):
    def __init__(self):
        Toplevel.__init__(self)
        self.title("Measurement model determination")

        global LtnNames
        global ObsToLtnTemp
        global isModelDone
        if not isModelDone:
            LtnNames = []
            ObsToLtnTemp = NumOfObsVarTemp * [-1]
        print(isModelDone)
        #  Input latent variable names
        label_frame1 = ttk.LabelFrame(self, text = "STEP1: Enter latent variables names")
        label_frame1.grid(row=0, column=0, columnspan=2, padx=10, pady=10,
                          sticky=(tk.N, tk.S, tk.E, tk.W))
        # label1 = tk.Label(label_frame1, text="Enter new latent variable names")
        # label1.grid(row=0, column=0, padx=10, pady=10)
        self.entry1 = tk.Entry(label_frame1)
        self.entry1.pack(fill="x", expand=True, side="top", padx=10, pady=5)
        self.entry1.bind("<Return>", self.add_ltn_name)
        #self.entry1.grid(row=0, column=0, padx=80, pady=10)
        button1 = tk.Button(label_frame1, text="Add latent variable", command=self.add_ltn_name)
        button1.pack(fill="x", expand=True, side="bottom", padx=10, pady=5)
        #button1.grid(row=1, column=0, padx=80, pady=10)

        #  link observed variables and latent variables
        label_frame2 = ttk.LabelFrame(self, text = "STEP2: Connect latent variables to observed variables")
        label_frame2.grid(row=0, column=2, columnspan=2,padx=10, pady=10,
                          sticky=(tk.N, tk.S, tk.E, tk.W))

        label1 = tk.Label(label_frame2, text="Select a latent variable")
        label1.grid(row=0, column=0, padx=10, pady=10)
        label2 = tk.Label(label_frame2, text="Select an observed variable")
        label2.grid(row=0, column=2, padx=10, pady=10)

        self.ltnname = tk.StringVar()
        self.ltn_combo = ttk.Combobox(label_frame2, width=15, textvariable=self.ltnname,
                                state="readonly")
        self.ltn_combo.grid(row=1, column=0, padx=10, pady=10)

        label_frame3 = ttk.LabelFrame(label_frame2, text="Connect")
        label_frame3.grid(row=0, column=1, rowspan=2, padx=10, pady=10,
                          sticky=(tk.N, tk.S, tk.E, tk.W))
        label3 = ttk.Label(label_frame3, text="Lavaan style: =~")
        label4 = ttk.Label(label_frame3, text="Mplus style: BY")
        label3.grid(column=0, row=0)
        label4.grid(column=0, row=1)
        label3.configure(justify="center")
        label4.configure(justify="center")

        button2 = tk.Button(label_frame3, text="OK", command=self.connect)
        button2.grid(row=2, column=0, padx=10, sticky="nsew")
        self.obsname = tk.StringVar()
        self.obs_combo = ttk.Combobox(label_frame2, width=15,
                                      textvariable=self.obsname,
                                      state="readonly")
        self.obs_combo["values"] = ObsNamesTemp
        self.obs_combo.grid(row=1, column=2, padx=10, pady=10)
        
        # Express user's input in Lavaan and mplus styles
        label5 = tk.Label(self, text="Lavaan view")
        label5.grid(row=1, column=0, padx=10)
        label6 = tk.Label(self, text="Mplus view")
        label6.grid(row=1, column=2, padx=10)
        
        self.txt1 = tk.Text(self, width=60, wrap=tk.WORD)
        self.txt1.grid(row=2, column=0)
        self.txt2 = tk.Text(self, width=60, wrap=tk.WORD)
        self.txt2.grid(row=2, column=2)
        scrl1=tk.Scrollbar(self,command=self.txt1.yview, width=16)
        scrl2=tk.Scrollbar(self,command=self.txt2.yview, width=16)
        scrl1.grid(row=2, column=1, sticky="nsew")
        scrl2.grid(row=2, column=3, sticky="nsw")
        self.txt1["yscrollcommand"] = scrl1.set
        self.txt2["yscrollcommand"] = scrl2.set

        label_frame4 = ttk.LabelFrame(self, text="STEP3: Reset or confirm")
        label_frame4.grid(row=3, column=2, columnspan=2, padx=10, pady=10,
                          sticky=(tk.N, tk.S, tk.E, tk.W))
        button3 = tk.Button(label_frame4, text="Reset",
                            height=1, width =10, command=self.reset)
        button3.pack(fill="both", expand=True, side="left", padx=10, pady=5)
        #button3.grid(row=0, column=1, padx=50, sticky=(tk.N, tk.S, tk.E, tk.W))
        button4 = tk.Button(label_frame4, text="OK", 
                            height=1, width =10, command=self.ok)
        button4.pack(fill="both", expand=True, side="right", padx=10, pady=5)
        #button4.grid(row=0, column=4, padx=50, sticky=(tk.N, tk.S, tk.E, tk.W))
        # -------------------------------------------------------
        # Displays existing models, if any.
        #-------------------------------------------------------        
        if LtnNames != []:
            self.ltn_combo["values"] = LtnNames
            self.obs_combo["values"] = ObsNamesTemp
            self.display()

    def add_ltn_name(self, event=None):
        global LtnNames
        global NumOfLtnVar
        EnteredName = self.entry1.get()
        LtnNames.append(EnteredName)
        NumOfLtnVar = len(LtnNames)
        self.ltn_combo["values"] = LtnNames
        self.display()
 
    def connect(self):
        global ObsToLtnTemp
        CurrentLtn = LtnNames.index(self.ltn_combo.get())
        CurrentObs = ObsNamesTemp.index(self.obs_combo.get())
        # Connects the observed variable and the latent variable
        ObsToLtnTemp[CurrentObs] = CurrentLtn
        self.display()

    def display(self):
        #  to display the connection status on the screen
        # in the style of mplus and Lavaan.
        self.txt1.delete("1.0", "end")
        self.txt2.delete("1.0", "end")
        global LtnNames
        
        for i in range(NumOfLtnVar):
            isFirst = True
            for j in range(NumOfObsVarTemp):
                if ObsToLtnTemp[j] == i:
                    if isFirst:
                        self.txt1.insert("%d.%d" % (i+1, 0),
                                         LtnNames[i] + " =~ NA * ")
                        self.txt2.insert("%d.%d" % (i+1, 0),
                                         LtnNames[i] + " BY 1 * ")
                        #  I used the if statement to omit the + sign only
                        #for the first observed variable.
                        self.txt1.insert("%d.end" % (i+1),
                                         ObsNamesTemp[j] + " ")
                    else:
                        self.txt1.insert("%d.end" % (i+1),
                                         " + " + ObsNamesTemp[j] + " ")
                    self.txt2.insert("%d.end" % (i+1), ObsNamesTemp[j] + " ")
                    isFirst = False
            self.txt1.insert("%d.end" % (i+1), "\n")
            self.txt2.insert("%d.end" % (i+1), "\n")
            #self.txt2.insert("%d.end" % (i+1), ";")
        for i in range(NumOfLtnVar):
            self.txt1.insert("%d.%d" % (i + NumOfLtnVar + 1, 0), 
                             LtnNames[i] + " ~~ 1 * "+ LtnNames[i] + "\n")
            self.txt2.insert("%d.%d" % (i + NumOfLtnVar  + 1, 0),
                             LtnNames[i] + " @1" + "\n")

    def reset(self):
        global ObsToLtnTemp
        global LtnNames
        global NumOfLtnVar
        ObsToLtnTemp = NumOfObsVarTemp * [-1]
        LtnNames = []
        NumOfLtnVar = 0
        self.txt1.delete("1.0", "end")
        self.txt2.delete("1.0", "end")

    def ok(self):
        # Ensure that all latent variables are specified.
        isAssigned = NumOfLtnVar * [False]
        for i in range(NumOfObsVarTemp):
            for j in range(NumOfLtnVar):
                if ObsToLtnTemp[i] == j:
                    isAssigned[j] = True

        if False in isAssigned:
            messagebox.showwarning("Unspecified latent variables",
              "Some latent variables are not connected to the observed variables.")

        else:
            global NumOfObsVar
            global ObsCovs
            global ObsNames
            global ObsToLtn
            global isModelDone
            NumOfObsVar = NumOfObsVarTemp
            ObsCovs = ObsCovsTemp[:]  
            ObsNames = ObsNamesTemp[:]
            ObsToLtn = ObsToLtnTemp[:]
            i = 0
            while i < len(ObsToLtn):
                if ObsToLtn[i] == -1:  # if i is not assigned
                    del ObsToLtn[i] # delete i 
                    del ObsNames[i] # delete i
                    ObsCovs = np.delete(ObsCovs, i, axis = 0)
                    ObsCovs = np.delete(ObsCovs, i, axis = 1)
                else:
                    i += 1
                    
            NumOfObsVar = len(ObsToLtn)
            
            isModelDone = True
            RunOption()
            self.destroy()


class RunOption(Toplevel):
    def __init__(self):
        Toplevel.__init__(self)
        self.title("Run options")
        
        # family-wise Confidence level
        frame1 = ttk.Labelframe(self,text="Confidence level")
        frame1.grid(row=0, column=0, padx=10, pady=5,
                    sticky=(tk.N, tk.S, tk.E, tk.W))
        lbl1 = ttk.Label(frame1, text="Choose the family-wise error rate or cumulative Type I error")
        lbl1.grid(row=0, column=0, columnspan=5, padx=10, pady=10, sticky=tk.W)
        self.radVar = tk.IntVar(value=2)  #  assigning a default value
        self.rad1 = tk.Radiobutton(frame1, text = "0.01", variable=self.radVar,
                                   value=1)
        self.rad1.grid(row=1, column=0)
        self.rad2 = tk.Radiobutton(frame1, text = "0.05", variable=self.radVar,
                                   value=2)
        self.rad2.grid(row=1, column=1)
        self.rad3 = tk.Radiobutton(frame1, text = "0.10", variable=self.radVar,
                                   value=3)
        self.rad3.grid(row=1, column=2)
        self.rad4 = tk.Radiobutton(frame1, text = "Other values",
                                   variable=self.radVar, value=4)
        self.rad4.grid(row=1, column=3)
        self.UserAlpha = tk.StringVar()
        self.Alpha_Chosen = ttk.Combobox(frame1, width=12,
                                         textvariable=self.UserAlpha,
                                         state="readonly")
        self.Alpha_Chosen["values"] = (0.001, 0.005, 0.02, 0.025, 0.03, 0.04,
                         0.06, 0.07, 0.08, 0.09, 0.11, 0.12, 0.13, 0.14, 0.15,
                         0.16, 0.17, 0.18, 0.19, 0.20)
        self.Alpha_Chosen.grid(row=1, column=4, padx=10, pady=10,)
        self.Alpha_Chosen.set("0.001") #  assigning a default value
        self.Alpha_Chosen.bind("<<ComboboxSelected>>",self.comborad1)
        
        # Cutoff point
        frame2 = ttk.Labelframe(self,text="Cutoff point")
        frame2.grid(row=1, column=0, padx=10, pady=5,
                    sticky=(tk.N, tk.S, tk.E, tk.W))
        lbl3 = tk.Label(frame2, text="For every factor pair, this program will test the statistical significance \n of the chi-square difference between the model with fixed factor correlations \n and the original model. Determine the value to be fixed.")
        lbl3.grid(row=0, column=0, columnspan=6, padx=10, pady=10, sticky=tk.W)
        lbl3.configure(justify="left")
        
        self.cut = tk.IntVar(value=2)  #  assigning a default value
        self.rad5 = tk.Radiobutton(frame2, text = "0.85", variable=self.cut, value=1)
        self.rad5.grid(row=1, column=0)
        self.rad6 = tk.Radiobutton(frame2, text = "0.90", variable=self.cut, value=2)
        self.rad6.grid(row=1, column=1)
        self.rad7 = tk.Radiobutton(frame2, text = "0.95", variable=self.cut, value=3)
        self.rad7.grid(row=1, column=2)
        self.rad8 = tk.Radiobutton(frame2, text = "1.0", variable=self.cut, value=4)
        self.rad8.grid(row=1, column=3)
        self.rad9 = tk.Radiobutton(frame2, text = "Other values", variable=self.cut,
                              value=5)
        self.rad9.grid(row=1, column=4)
        self.UserCut = tk.StringVar()
        self.Cut_Chosen = ttk.Combobox(frame2, width=12, textvariable=self.UserCut,
                                  state="readonly")
        self.Cut_Chosen["values"] = (0.80, 0.81, 0.82, 0.83, 0.84, 0.86, 0.87, 0.88,
                  0.89, 0.91, 0.92, 0.93, 0.94, 0.96, 0.97, 0.98, 0.99)
        self.Cut_Chosen.grid(row=1, column=5, padx=10, pady=10, sticky=tk.W)       
        self.Cut_Chosen.set("0.80") #  assigning a default value
        self.Cut_Chosen.bind("<<ComboboxSelected>>",self.comborad2)

        # Decimal place
        frame3 = ttk.Labelframe(self,text="Decimal places")
        frame3.grid(row=2, column=0, padx=10, pady=5,
                    sticky=(tk.N, tk.S, tk.E, tk.W))
        lbl4 = tk.Label(frame3, text="How many decimal places should be used in the output? (e.g., 3.456789)")
        lbl4.grid(row=0, column=0, columnspan=6, padx=10, pady=10, sticky=tk.W)
        lbl3.configure(justify="left")
        
        self.dec = tk.IntVar(value=2)
        rad10 = tk.Radiobutton(frame3, text = "3.46", variable=self.dec, value=1)
        rad10.grid(row=1, column=0)
        rad11 = tk.Radiobutton(frame3, text = "3.457", variable=self.dec,
                               value=2)
        rad11.grid(row=1, column=1)
        rad12 = tk.Radiobutton(frame3, text = "3.4568", variable=self.dec,
                               value=3)
        rad12.grid(row=1, column=2)
        rad13 = tk.Radiobutton(frame3, text = "3.45679", variable=self.dec,
                               value=4)
        rad13.grid(row=1, column=3)
        rad14 = tk.Radiobutton(frame3, text = "3.456789", variable=self.dec,
                              value=5)
        rad14.grid(row=1, column=4)
        
        # Comments
        frame4 = ttk.Labelframe(self,text="Comments")
        frame4.grid(row=3, column=0, padx=10, pady=5,
                    sticky=(tk.N, tk.S, tk.E, tk.W))
        lbl5 = tk.Label(frame4, text="Do you want to include methodological comments in the output?")
        lbl5.grid(row=0, column=0, columnspan=6, padx=10, pady=10, sticky=tk.W)
        lbl5.configure(justify="left")
        
        self.com = tk.IntVar(value=1)
        rad15 = tk.Radiobutton(frame4, text = "Yes", variable=self.com, value=1)
        rad15.grid(row=1, column=0, padx = 10)
        rad16 = tk.Radiobutton(frame4, text = "No", variable=self.com, value=2)
        rad16.grid(row=1, column=1, padx = 10)

        # Ok  button
        btn1 = tk.Button(self, text="OK", height=1, width =10, command=self.ok)
        btn1.grid(row=5, column=0, padx=10, pady=10)
        
        self.iconbitmap('taichi.ico')

    def comborad1(self, event):
        self.rad1.deselect()
        self.rad2.deselect()
        self.rad3.deselect()
        self.rad4.select()

    def comborad2(self, event):
        self.rad5.deselect()
        self.rad6.deselect()
        self.rad7.deselect()
        self.rad8.deselect()
        self.rad9.select()
        
    def ok(self):
        #  Process the first module to change the global variable Alpha.
        global Alpha
        alp = self.radVar.get()
        if alp == 1:
            Alpha = 0.01
        elif alp == 2:
            Alpha = 0.05
        elif alp == 3:
            Alpha = 0.10
        else:
            Alpha = float(self.UserAlpha.get())
        
        #  Process the second module to change the global variable Cutoff.
        global Cutoff
        ctf = self.cut.get()
        if ctf == 1:
            Cutoff = 0.85
        elif ctf == 2:
            Cutoff = 0.90
        elif ctf == 3:
            Cutoff = 0.95
        elif ctf == 4:
            Cutoff = 1.0
        else:
            Cutoff = float(self.UserCut.get())
        
        #  Process the third module to change the global variable Dcm.
        global Dcm
        decimal = self.dec.get()
        if decimal == 1:
            Dcm = 2
        elif decimal == 2:
            Dcm = 3
        elif decimal == 3:
            Dcm = 4
        elif decimal == 4:
            Dcm = 5
        else:
            Dcm = 6
        
        #  Process the module to change the global variable isCmt.
        global isCmt
        comment = self.com.get()
        if comment == 1:
            isCmt = True
        else:
            isCmt = False

        
        Output()
        
        # close the window
        self.destroy()

class Output(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.title("Output")
        TWIDTH = 250
        THEIGHT = 80
        #creating a tab
        tabControl = ttk.Notebook(self)
        tab1 = ttk.Frame(tabControl)
        tab2 = ttk.Frame(tabControl)
        tab3 = ttk.Frame(tabControl)
        tab4 = ttk.Frame(tabControl)
        tab5 = ttk.Frame(tabControl)
        tabControl.add(tab1, text = "  Model fit  ")
        tabControl.pack(expand=1, fill="both", ipadx=20)
        tabControl.add(tab2, text = "  Parameter estimates  ")
        tabControl.pack(expand=1, fill="both")
        tabControl.add(tab3, text = "  Discriminant validity  ")
        tabControl.pack(expand=1, fill="both")
        tabControl.add(tab4, text = "  Reliability estimates  ")
        tabControl.pack(expand=1, fill="both")
        tabControl.add(tab5, text = "  Correlation table  ")
        tabControl.pack(expand=1, fill="both")
        
        
        # -------------------------------------------------------
        # Perform various calculations
        #-------------------------------------------------------
        global NumOfParam
        global NumOfComparison
        global FixPhi
        global InitialCutoff
        # NumOfParam represents the number of unknown parameters.
        # Because of the division by 2, I added the int () otherwise the computer
        # recognized NumOfParam as float instead of int.
        NumOfParam = int(2 * NumOfObsVar + (NumOfLtnVar) * (NumOfLtnVar - 1) / 2)
        NumOfComparison = int(NumOfLtnVar * (NumOfLtnVar - 1) / 2)
        InitialCutoff = Cutoff
        # Converts observed data into observed covariances.
        # Initial guess is that all parameters will be 1.0 (float)
        InitialGuess = np.ones(NumOfParam)
        FixPhi = -1  # -1 means that no phi is fixed
        # The core function used to estimate mle is scipy.optimize.minimize
        # I have confirmed the results of several methods of the minimize function.
        # Of the methods that do not require jacobian, the default option took
        # the fastest time. Perhaps the default option is BFGS, not 100% sure.
        Result = minimize(ml, InitialGuess)
        Stats = get_stats(Result.x, Result.fun, NumOfParam)
        Params = x_to_param(Result.x)
        Lambdas = Params["Lambdas"]
        Deltas = Params["Deltas"]
        PhiMatrix = Params["PhiMatrix"] 
        StdLambdas = get_stdlambdas(Lambdas, Deltas)
        SS = get_sscor()
        Rel = get_reliability(Lambdas, Deltas, PhiMatrix)
        
        # Convert the x values obtained as a result of the minimize function
        # back to the parameters we use, and print them.
        
        # -------------------------------------------------------
        # print Model fit to the screen
        #-------------------------------------------------------
        self.txt1 = tk.Text(tab1, width=TWIDTH, height=THEIGHT, wrap=tk.WORD)
        self.txt1.pack(expand=1, fill="both")
        self.txt1.insert("1.0", "Model fit summary:\n")
        if Stats["CFI"] > .95 and Stats["TLI"] > .95 and Stats["RMSEA"] < .06 \
        and Stats["SRMR"] < .08:
            self.txt1.insert("end", "The fit indexes of the measurement model satisfied the cutoff \
criteria of Hu and Bentler (1999) with CFI " + format(Stats["CFI"], ("4.%df")%(Dcm))
              + ", TLI " + format(Stats["TLI"], ("4.%df")%(Dcm))
              + ", RMSEA", format(Stats["RMSEA"], ("4.%df")%(Dcm))
              + " and SRMR" + format(Stats["SRMR"], ("4.%df")%(Dcm)) + ".\n")
        elif Stats["CFI"] > .95 or Stats["TLI"] > .95 or Stats["RMSEA"] < .06 \
            or Stats["SRMR"] < .08:
            self.txt1.insert("end", "The fit indexes of the measurement model satisfied the cutoff \
criteria of Hu and Bentler (1999) with\n")
            if Stats["CFI"] > .95:
                self.txt1.insert("end"," CFI (" + format(Stats["CFI"], ("4.%df")%(Dcm)) + ")")
            if Stats["TLI"] > .95:
                self.txt1.insert("end"," TLI (" + format(Stats["TLI"], ("4.%df")%(Dcm)) + ")")
            if Stats["RMSEA"] < .06:
                self.txt1.insert("end", " RMSEA ("+ format(Stats["CFI"], ("4.%df")%(Dcm)) + ")")
            if Stats["SRMR"] < .08:
                self.txt1.insert("end", " SRMR (" + format(Stats["SRMR"], ("4.%df")%(Dcm)) + ")")
            self.txt1.insert("end",", but ")
            if Stats["CFI"] < .95:
                self.txt1.insert("end", " CFI (" + format(Stats["CFI"], ("4.%df")%(Dcm)) + ")")
            if Stats["TLI"] < .95:
                self.txt1.insert("end", " TLI (" + format(Stats["TLI"], ("4.%df")%(Dcm)) + ")")
            if Stats["RMSEA"] > .06:
                self.txt1.insert("end", " RMSEA ("+ format(Stats["RMSEA"], ("4.%df")%(Dcm)) + ")")
            if Stats["SRMR"] > .08:
                self.txt1.insert("end", " SRMR ("+ format(Stats["SRMR"], ("4.%df")%(Dcm)) + ")")
            self.txt1.insert("end"," did not.\n")
        else:
            self.txt1.insert("end","The fit indexes of the measurement model did not satisfy the\
cutoff criteria of Hu and Bentler (1999) with\n CFI " + format(Stats["CFI"], ("4.%df")%(Dcm))
                  + ", TLI " + format(Stats["TLI"], ("4.%df")%(Dcm))
                  + ", RMSEA " + format(Stats["RMSEA"], ("4.%df")%(Dcm))
                  + " and SRMR " + format(Stats["SRMR"], ("4.%df")%(Dcm)) + ".\n")
        
        self.printline(self.txt1)
        self.table(self.txt1, "Model", "Number of parameters", "Chi-square", 
                   "degree of freedom(Df)", "p", "Chi-square/Df", 
                   "Current Model ",
                   "%d"%(NumOfParam), format(Stats["Chi2"], ("4.%df")%(Dcm)), 
                   "%d"%(Stats["Df"]), format(Stats["PValue"], ("4.%df")%(Dcm)),
                   format(Stats["Chi2DivDf"], ("4.%df")%(Dcm)), 
                   "Null Model", "%d"%(NumOfObsVar), 
                   format(Stats["Chi2Null"], ("4.%df")%(Dcm)),
                   "%d"%(Stats["DfNull"]), 
                   format(Stats["PValueNull"], ("4.%df")%(Dcm)), 
                   format(Stats["Chi2DivDfNull"], ("4.%df")%(Dcm)),
                   isList=False, row=3)
        
        self.txt1.insert("end", "\n" + "Model fit indexes" + "\n")
        self.table(self.txt1, "NFI (Normed Fit Index or Bentler-Bonnet Index): ", 
                   format(Stats["NFI"], ("4.%df")%(Dcm)),
                   "TLI (Tucker-Lewis Index or Non-normed Fit Index): ",
                   format(Stats["TLI"], ("4.%df")%(Dcm)),
                   "CFI (Comparative Fit Index): ",
                   format(Stats["CFI"], ("4.%df")%(Dcm)),
                   "IFI (Incremental Fit Index): ",
                   format(Stats["IFI"], ("4.%df")%(Dcm)),
                   "RMSEA (Root Mean Squared Error of Approximation): ",
                   format(Stats["RMSEA"], ("4.%df")%(Dcm)),
                   "SRMR (Standardized Root Mean Squared Residual): ",
                   format(Stats["SRMR"], ("4.%df")%(Dcm)),
                   "AIC (Akaike Information Criterion): ",
                   format(Stats["AIC"], ("4.%df")%(Dcm)),
                   "BIC (Bayesian Information Criterion): ",
                   format(Stats["BIC"], ("4.%df")%(Dcm)),
                   "SABIC (Sample-size Adjusted BIC): ",
                   format(Stats["SABIC"], ("4.%df")%(Dcm)), isList=False, row=9)

        self.printline(self.txt1)
        
        # -------------------------------------------------------
        # print parameter estimates to the screen
        #-------------------------------------------------------
        self.txt2 = tk.Text(tab2, width=TWIDTH, height=THEIGHT, wrap=tk.WORD)
        self.txt2.pack(expand=1, fill="both")
        self.txt2.insert("1.0", "Parameter estimates:\n")

        self.txt2.insert("end", "Regression coefficients between observed and latent variables (lambda in the LISREL notation)" + "\n")
        
        # Work to show the coefficients in a table format
        global Coeffs
        Coeffs = [[""]*3 for i in range(NumOfObsVar + 1)]
        Coeffs[0][0] = " "
        Coeffs[0][1] = "Coeffcient(lambda)"
        Coeffs[0][2] = "Standardized Coefficient"
        for i in range(NumOfObsVar):
            Coeffs[i+1][0] = ObsNames[i] + "<--" + LtnNames[ObsToLtn[i]] + ": "
            Coeffs[i+1][1] = format(Lambdas[i], ("4.%df")%(Dcm))
            Coeffs[i+1][2] = format(StdLambdas[i], ("4.%df")%(Dcm))
        Coeff1 = sum(Coeffs, [])  # making a one-dimensional list
        self.table(self.txt2, Coeff1, isList=True, row=NumOfObsVar + 1)
        self.printline(self.txt2)
        
        # Print Error variances
        self.txt2.insert("end", "Error variances of observed variables (delta in the LISREL notation)" + "\n")
        global Resids
        Resids = [[""]*2 for i in range(NumOfObsVar)]
        for i in range(NumOfObsVar):
            Resids[i][0] = ObsNames[i] + "  : "
            Resids[i][1] = format(Deltas[i], ("4.%df")%(Dcm))
        Resid1 = sum(Resids, [])  # making a one-dimensional list
        self.table(self.txt2, Resid1, isList=True, row=NumOfObsVar)
        self.printline(self.txt2)
        
        # Print covariances
        self.txt2.insert("end", "Covariances/correlations between latent variables (phi in the LISREL notation)" + "\n")
        global Cors
        Cors = [[""]*2 for i in range(NumOfComparison)]
        count = 0
        for i in range(NumOfLtnVar):
            for j in range(i + 1, NumOfLtnVar):
                Cors[count][0] = LtnNames[i] + " <-> " + LtnNames[j]
                Cors[count][1] = format(PhiMatrix[i][j], ("4.%df")%(Dcm))
                count += 1
        Cor1 =  sum(Cors, [])  # making a one-dimensional list
        self.table(self.txt2, Cor1, isList=True, row=NumOfComparison)
        self.printline(self.txt2)       
        # -------------------------------------------------------
        # print discriminant validity to the screen
        #-------------------------------------------------------
        
        # -------------------------------------------------------
        # print discriminant validity summary
        #-------------------------------------------------------
        self.txt3 = tk.Text(tab3, width=TWIDTH, height=THEIGHT, wrap=tk.WORD)
        self.txt3.pack(expand=1, fill="both")
        AllPairs = chi2diff_allpairs(Result.x, Result.fun, NumOfParam, Cutoff)
        CompCutoff = AllPairs["CompCutoff"]
        self.txt3.insert("1.0", "Discriminant valdity summary:\n")
        if "Greater" in CompCutoff or "NotSig" in CompCutoff:
            self.txt3.insert("end", "There is a problem with the discrimination validity.\n")
            if "Greater" in CompCutoff:
                self.txt3.insert("end", "Some of the correlations were larger than the cutoff " +
                                 format(Cutoff,".2f") +" you set.\n")
                count = 0
                for i in range(NumOfLtnVar):
                    for j in range(i + 1, NumOfLtnVar):
                        if CompCutoff[count] == "Greater":
                            self.txt3.insert("end", LtnNames[i] + " <-> " 
                                             + LtnNames[j] + " (" +
                                             format(PhiMatrix[i][j], ("4.%df")%(Dcm)) + ")\n")
                        count += 1
            if "NotSig" in CompCutoff:
                self.txt3.insert("end", "Some of the correlations were not significantly different from the cutoff "
                                 + format(Cutoff, ".2f") + " you set.\n")
                count = 0
                for i in range(NumOfLtnVar):
                    for j in range(i + 1, NumOfLtnVar):
                        if CompCutoff[count] == "NotSig":
                            self.txt3.insert("end", LtnNames[i] + " <-> " +
                                             LtnNames[j] + " (" + 
                                             format(PhiMatrix[i][j], ("4.%df")%(Dcm)) + ")\n")
                        count += 1
        else:
            self.txt3.insert("end", "No problems were found in the discriminant validity. All \
correlations were significantly smaller than the cutoff you set.\n")
        self.printline(self.txt3)
        global DVs
        DVs = [[""]*5 for i in range(NumOfComparison + 2)]
        DVs[0][0] = "Fixed pair"
        DVs[0][1] = "Correlation before fixed"
        DVs[0][2] = "Fixed correlation value"
        DVs[0][3] = "Chi-square"
        DVs[0][4] = "CFI"
        DVs[NumOfComparison + 1][0] = "(Original model)"
        DVs[NumOfComparison + 1][1] = DVs[NumOfComparison + 1][2] = "-"
        DVs[NumOfComparison + 1][3] = format(Stats["Chi2"], ("4.%df")%(Dcm))
        DVs[NumOfComparison + 1][4] = format(Stats["CFI"], ("4.%df")%(Dcm))
        count = 1  # 
        for i in range(NumOfLtnVar):
            for j in range(i + 1, NumOfLtnVar):
                DVs[count][0] = LtnNames[i] + " <-> " + LtnNames[j]
                DVs[count][1] = format(PhiMatrix[i][j], ("4.%df")%(Dcm))
                if PhiMatrix[i][j] > 0:
                    DVs[count][2] = format(Cutoff, "4.2f")
                else:
                    DVs[count][2] = format(Cutoff * (-1), "4.2f")
                DVs[count][3] = format(AllPairs["Chi2Fix"][count - 2], ("4.%df")%(Dcm))
                DVs[count][4] = format(AllPairs["CFIFix"][count - 2], ("4.%df")%(Dcm))
                count += 1
        DV1 =  sum(DVs, [])  # making a one-dimensional list
        self.table(self.txt3, DV1, isList=True, row=NumOfComparison + 2)
        
        if isCmt:
            Cmt_DV = "\nNumerous techniques for discriminant validity have been \
used, but there has been insufficient discussion of which technique is best. \
The results presented above represent an exemplary procedure based on our \
study. It has been known for a long time that the technique of testing the \
chi - square difference between the model with the correlation between the \
two factors fixed at 1.0 and the original model. Our new suggestion is that \
fixing the correlation to 1.0 is not enough and you should use a cutoff point \
lower than 1.0. Running proposed pair-wise model comparisons in conventional SEM \
software is time-consuming and error-prone. That's why we created this program. In the \
following section, we present some logical problems with the techniques that \
have been used to evaluate discriminant validity. Although correlations are \
key to discriminant validity, the most commonly used techniques are only \
indirectly related to correlations. See the following article for\
 details.\n  Discriminant validity: what is it and how to assess it, \
working paper\n"
        self.get_cmt(self.txt3, Cmt_DV)    
        self.printline(self.txt3)   
        
        # -------------------------------------------------------
        # print discriminant validity : AVE
        #-------------------------------------------------------
        # Calculating AVE
        NumSum = NumOfLtnVar * [0]
        DenSum = NumOfLtnVar * [0]
        AVE = NumOfLtnVar * [0]
        for i in range(NumOfObsVar):
            for j in range(NumOfLtnVar):
                NumSum[ObsToLtn[i]] += Lambdas[i]**2
                DenSum[ObsToLtn[i]] += Lambdas[i]**2 + Deltas[i]
        for j in range(NumOfLtnVar):
            AVE[j] = NumSum[j] / DenSum[j]
            
        # Displaying AVE
        self.txt3.insert("end", "Commonly used but inapprorpriate techniques: AVE/SV\n")
        global AVEs
        AVEs = [[""]*5 for i in range(NumOfComparison + 1)]
        AVEs[0][0] = " "
        AVEs[0][1] = "Shared variance (SV)"
        AVEs[0][2] = "AVE1"
        AVEs[0][3] = "AVE2"
        AVEs[0][4] = "SV < max(AVE1, AVE2)?"
        count = 1
        for i in range(NumOfLtnVar):
            for j in range(i + 1, NumOfLtnVar):
                AVEs[count][0] = LtnNames[i] + "<->" + LtnNames[j]
                AVEs[count][1] = format(PhiMatrix[i][j]**2, ("4.%df")%(Dcm))
                AVEs[count][2] = format(AVE[i], ("4.%df")%(Dcm))
                AVEs[count][3] = format(AVE[j], ("4.%df")%(Dcm))
                if PhiMatrix[i][j]**2 < max(AVE[i], AVE[j]):
                    AVEs[count][4] = "Pass" 
                else:
                    AVEs[count][4] = "Fail"
                count += 1
        AVE1 =  sum(AVEs, [])  # making a one-dimensional list
        self.table(self.txt3, AVE1, isList=True, row=NumOfComparison + 1)
                
        if isCmt:
            Cmt_AVE = "\nAVE / SV is the criterion that the maximum value of \
the Average Variance Extracted (AVE) of the two latent variables must be \
greater than the shared variance (SV), the square of the factor correlation \
(CFA) between two latent variables, . This technique is the most commonly used technique \
for evaluating discriminant validity in many fields including marketing. \
However, AVE is UNRELATED with discriminant validity. It is instead an \
indicator of item reliability. Most studies are applying this criterion \
differently from the original proposal of Fornell and Larcker (1984). Meeting \
this criterion is very difficult if it is applied 'correctly'. An irony is \
that its high false positive rate has been interpreted positively as a sign \
of the statistical power, and may have contributed to the popularity of \
this technique.\n"
            self.get_cmt(self.txt3, Cmt_AVE)    
        self.printline(self.txt3)   
        # -------------------------------------------------------
        # print discriminant validity : HTMT
        #-------------------------------------------------------
        # Calculating disattenuated correlations
        CFRel = Rel["CFRel"]
        SubRel = Rel["SubRel"]
        ConRel = Rel["ConRel"]
        TauEqRel = Rel["TauEqRel"]
        ParRel = Rel["ParRel"]
        SSCor = SS["SSCor"]
        # Disattenuated correlation using congeneric reliability      
        DCCR = NumOfComparison * [0]
        # Disattenuated correlation using tau-equivalent reliability
        DCTR = NumOfComparison * [0]
        # # Disattenuated correlation using parallel reliability
        DCPR = NumOfComparison * [0]
        # Obtaining disattneuated correlations
        count = 0
        for i in range(NumOfLtnVar):
            for j in range(i + 1, NumOfLtnVar):
                DCCR[count] = SSCor[i][j] / sqrt(ConRel[i] * ConRel[j])
                DCTR[count] = SSCor[i][j] / sqrt(TauEqRel[i] * TauEqRel[j])
                DCPR[count] = SSCor[i][j] / sqrt(ParRel[i] * ParRel[j])
                count += 1
        # Displaying
        self.txt3.insert("end", "Commonly used but inapprorpriate techniques: HTMT\n")
        self.txt3.insert("end", "Error adjusted correlation estimates\n")
        global DCs
        DCs = [[""]*5 for i in range(NumOfComparison + 1)]
        DCs[0][0] = " "
        DCs[0][1] = "CFA"
        DCs[0][2] = "DCCR"
        DCs[0][3] = "DCTR"
        DCs[0][4] = "DCPR"
        count = 1
        for i in range(NumOfLtnVar):
            for j in range(i + 1, NumOfLtnVar):
                DCs[count][0] = LtnNames[i] + "<->" + LtnNames[j]
                DCs[count][1] = format(PhiMatrix[i][j], ("4.%df")%(Dcm))
                DCs[count][2] = format(DCCR[count - 1], ("4.%df")%(Dcm))
                DCs[count][3] = format(DCTR[count - 1], ("4.%df")%(Dcm))
                DCs[count][4] = format(DCPR[count - 1], ("4.%df")%(Dcm))
                count += 1
        DC1 =  sum(DCs, [])  # making a one-dimensional list
        self.table(self.txt3, DC1, isList=True, row=NumOfComparison + 1)
                
        if isCmt:
            Cmt_HTMT = "\nHTMT: Hetero-trait Mono-trait, DCCR: Disattenuated correlation \
using congeneric reliability, DCTR: Disattenuated correlation using \
tau-equivalent reliability, DCPR: Disattneuated correlation using parallel\
reliability. Recently introduced as a new technique, HTMT is actually a very \
old technique called DCPR and it is the least accurate method among several \
disattenuated correlations. Disattenuated correlation techniques have some \
disadvantages compared to CFA.\n"
            self.get_cmt(self.txt3, Cmt_HTMT)
        self.printline(self.txt3)   
        # -------------------------------------------------------
        # print discriminant validity : Cross-loadings
        #-------------------------------------------------------
        # Obtaining factor Pattern and structure matrix
        PatnMatrix = get_patnmatrix(get_stdlambdas(Lambdas, Deltas))
        StrucMatrix = np.matmul(PatnMatrix, PhiMatrix)
        MinPatns = NumOfLtnVar * [1]
        # Find the minimum value of loading belonging to latent variable j
        for i in range(NumOfObsVar):
            for j in range(NumOfLtnVar):
                if ObsToLtn[i] == j and abs(PatnMatrix[i][j]) < abs(MinPatns[j]):
                    MinPatns[j] = PatnMatrix[i][j]
        # Displaying factor pattern matrix
        self.txt3.insert("end", "Commonly used but inapprorpriate techniques: Cross-loadings\n")
        self.txt3.insert("end", "Factor Pattern Coefficients\n")
        global FPs
        FPs = [[""]*(NumOfLtnVar + 1) for i in range(NumOfObsVar + 1)]
        FPs[0][0] = " "
        for i in range(NumOfLtnVar):
            FPs[0][i + 1] = LtnNames[i]
        for i in range(NumOfObsVar):
            FPs[i + 1][0] = ObsNames[i]
            for j in range(NumOfLtnVar):
                FPs[i + 1][j + 1] = format(PatnMatrix[i][j], ("4.%df")%(Dcm))
        FP1 = sum(FPs, [])
        self.table(self.txt3, FP1, isList=True, row=NumOfObsVar + 1)
        self.txt3.insert("end", "\n")
        # Displaying factor structure matrix
        self.txt3.insert("end", "Factor Structure Coefficients\n")
        global FSs
        FSs = [[""]*(NumOfLtnVar + 1) for i in range(NumOfObsVar + 1)]
        FSs[0][0] = " "
        for i in range(NumOfLtnVar):
            FSs[0][i + 1] = LtnNames[i]
        for i in range(NumOfObsVar):
            FSs[i + 1][0] = ObsNames[i]
            for j in range(NumOfLtnVar):
                FSs[i + 1][j + 1] = format(StrucMatrix[i][j], ("4.%df")%(Dcm))
        FS1 = sum(FSs, [])
        self.table(self.txt3, FS1, isList=True, row=NumOfObsVar + 1)
        # Displaying cross-loadings
        self.txt3.insert("end", "Cross-loadings according to row comparisons\n")
        Cmt_rc = "(The absolute value of the loading between a factor and an item \
belonging to that factor must be greater than those of the loadings between \
the item and any other factors)\n"
        self.get_cmt(self.txt3, Cmt_rc)
        self.txt3.insert("end", "[None]\n")
        self.txt3.insert("end", "\nCross-loadings according to column comparisons\n")
        Cmt_cc = "(The absolute value of the loading between a factor and an item \
belonging to that factor must be greater than those of the loadings between \
the factor and any items not belonging to that factor.)\n"
        self.get_cmt(self.txt3, Cmt_cc)

        
        for j in range(NumOfLtnVar):
            for i in range(NumOfObsVar):
                if ObsToLtn[i] != j and \
                abs(StrucMatrix[i][j]) > abs(MinPatns[j]):
                    self.txt3.insert("end", "The correlation between " +
                                     ObsNames[i] + " and " + LtnNames[j] + " " +
                                     format(StrucMatrix[i][j], ("4.%df")%(Dcm)) +
                                     " is greater than " + 
                                     format(MinPatns[j],  ("4.%df")%(Dcm)) + "\n")
        if isCmt:
            Cmt_cross = "\nFirst, cross-loadings have historically been a separate concept \
from discriminant validity. Second, there is no cross-loading in the sense \
that we commonly refer to. Loadings refer to either pattern coefficients \
(i.e., regression weights) and structure coefficients (i.e., correlation), \
depending on the context. Most SEM software does not show structure \
coefficients, and loadings in CFA almost always refer to pattern coefficients.\
 However, cross-loadings do no exist in pattern coefficients. It is unclear\
 what \
 criteria must be met to be classified as cross-loading, and any existing \
 criteria have a logical problem. Cross-loadings by row comparisons do not \
 exist in a typical CFA. On the contrary, column comparisons \
 classify too many loadings as cross-loadings.\n"
            self.get_cmt(self.txt3, Cmt_cross)
        self.printline(self.txt3)   
        # -------------------------------------------------------
        # print reliability estimates to the screen
        #-------------------------------------------------------
        self.txt4 = tk.Text(tab4, width=TWIDTH, height=THEIGHT, wrap=tk.WORD)
        self.txt4.pack(expand=1, fill="both")
        self.txt4.insert("1.0", "Reliability estimates:\n")
        self.txt4.insert("end", "Correlated factors reliability of the model is " +
                         format(CFRel,  ("4.%df")%(Dcm)) + ".\n" +
                         "Various reliability estimates for each latent variable are as follows.\n\n")
        global Rels
        Rels = [[""]*5 for i in range(NumOfLtnVar + 1)]
        Rels[0][0] = " "
        Rels[0][1] = "Subtest reliability"
        Rels[0][2] = "Congeneric reliability"
        Rels[0][3] = "Tau-equivalent reliability"
        Rels[0][4] = "Parallel reliability"
        count = 1
        for i in range(NumOfLtnVar):
            Rels[count][0] = LtnNames[i]
            Rels[count][1] = format(SubRel[i], ("4.%df")%(Dcm))
            Rels[count][2] = format(ConRel[i], ("4.%df")%(Dcm))
            Rels[count][3] = format(TauEqRel[i], ("4.%df")%(Dcm))
            Rels[count][4] = format(ParRel[i], ("4.%df")%(Dcm))
            count += 1
        Rel1 =  sum(Rels, [])  # making a one-dimensional list
        self.table(self.txt4, Rel1, isList=True, row=NumOfLtnVar + 1)
        if isCmt:
            Cmt_Rel = "Correlated factor reliability is a multidimensional reliability \
coefficient that is not recommended to be presented alone without subtest \
reliability coefficients. Subtest reliability and congeneric reliability have \
analogous formulas, both of which are often  referred to as composite \
reliability. The difference between the two is that  the former is calculated \
from the values obtained from the correlated factors  model and the latter is \
calculated from the values obtained from the unidimimensional model. \
Tau-equivalent reliability is often referred to as  the historically \
incorrect name Cronbach's alpha. Most experts disagree with its use because \
it tends to underestimate reliability and is less accurate. Parallel reliability\
 (i.e., standardized alpha) is less accurate than tau-equivalent reliability. See the \
following article for more information.\n  Cho, E. (2016). Making reliability \
reliable: A systematic approach\
 to reliability coefficients. Organizational Research Methods, 19(4), 651-682.\n"
            self.get_cmt(self.txt4, Cmt_Rel)
        self.printline(self.txt4)
        # -------------------------------------------------------
        # print correlation tables to the screen
        #-------------------------------------------------------
        self.txt5 = tk.Text(tab5, width=TWIDTH, height=THEIGHT, wrap=tk.WORD)
        self.txt5.pack(expand=1, fill="both")
        self.txt5.insert("1.0", "A copy-and-paste correlation table:\n")
        # Obtaining two correlation-related tables
        Matrix1 = np.zeros((NumOfLtnVar, NumOfLtnVar))
        PhiVector = Params["PhiVector"]
        # Matrix1: Scale score correlations are below diagonals,
        # Matrix2: P-values of scale score correlations are below diagonals,
        # factor correlations above diagonals, congeneric reliability on diagonals.
        count = 0
        for i in range(NumOfLtnVar):
            Matrix1[i][i] = SubRel[i]
            for j in range(i+1, NumOfLtnVar):
                Matrix1[i][j] = PhiVector[count]
                Matrix1[j][i] = SSCor[i][j]
                count += 1
        # Displyaing the first table
        self.txt5.insert("end", "Sub-diagonals: scale score correlations, \
diagonals: sub-test reliability, super-diagonals: factor correlation (CFA)\n")
        global Tbls
        Tbls = [[""]*(NumOfLtnVar + 1) for i in range(NumOfLtnVar + 1)]
        Tbls[0][0] = " "
        for i in range(NumOfLtnVar):
            Tbls[0][i + 1] = Tbls[i + 1][0] = LtnNames[i]
            for j in range(NumOfLtnVar):
                Tbls[i + 1][j + 1] = format(Matrix1[i][j], ("4.%df")%(Dcm))
        Tbl1 =  sum(Tbls, [])  # making a one-dimensional list
        self.table(self.txt5, Tbl1, isList=True, row=NumOfLtnVar + 1)
        if isCmt:
            Cmt_tbl = "In most behavior studies, a correlation matrix is presented, \
but the super-diagonals are hardly utilized. If factor correlation estimates \
are included at the top of the diagonal, \
more information can be provided without consuming additional space. The \
tables above are an example."
            self.get_cmt(self.txt5, Cmt_tbl)

        # Change icon
        self.iconbitmap('taichi.ico')
        
        global AllOutput
        AllOutput = self.txt1.get('1.0', 'end') + self.txt2.get('1.0', 'end') + self.txt3.get('1.0', 'end') + self.txt4.get('1.0', 'end') + self.txt5.get('1.0', 'end') 
        # -------------------------------------------------------
        # Save the state that the output is finished.
        #-------------------------------------------------------
        global isOutputDone
        isOutputDone = True
        
            

    def printline(self, where):
        LINE = 100
        where.insert("end", "\n")
        for i in range(LINE):
            where.insert("end", "-")
        where.insert("end", "\n")


    def table(self, where, *arg, isList=True, row=1):
        if isList:
            Cnt = arg[0]
        else:
            Cnt = arg
        col = round(len(Cnt)/row)
        colMax = col * [0]  # Length of the longest string in the columns
        for i in range(row):
            for j in range(col):
                if len(Cnt[i * col + j]) > colMax[j]:
                    colMax[j] = len(Cnt[i * col + j])
                    
        # The length of the column is the length of the longest string 
        # in the column plus two spaces.            
        for i in range(row):
            for j in range(col):
                where.insert("end", Cnt[i * col + j])
                nspace = colMax[j] - len(Cnt[i * col + j]) + 2
                for k in range(nspace):
                    where.insert("end", " ")
                
            where.insert("end", "\n")
        where.insert("end", "\n")
        
    def get_cmt(self, where, string):
        LINE = 100
        NumberOfLine = ceil(len(string)/LINE)
        NewTexts = NumberOfLine * [" "]
        current = 0
        for i in range(NumberOfLine):
            search = current + LINE
            if search < len(string):
                while string[search] != " ":
                    search -= 1
                if i > 0:
                    NewTexts[i] = NewTexts[i - 1] + string[current:search] + "\n"
                else: 
                    NewTexts[i] = string[current:search] + "\n"
                current = search
            else:
                if i > 0:
                    NewTexts[i] = NewTexts[i - 1] + string[current:len(string) - 1] + "\n"
                else:
                    NewTexts[i] = string[current:len(string) -1] + "\n"
        where.insert("end", NewTexts[NumberOfLine - 1])
'''
    def saveoutput(self):
        self.filename = tk.filedialog.asksaveasfilename(defaultextension='.txt',
                                                        filetypes = [('all files', '.*'), ('text files', '.txt')])
        f = open(self.filename, 'w')
        a = self.txt1.get('1.0', 'end') + self.txt2.get('1.0', 'end') + self.txt3.get('1.0', 'end') + self.txt4.get('1.0', 'end') + self.txt5.get('1.0', 'end') 
        f.write(a)
        f.close()
        tk.messagebox.showinfo('Save', 'Output Saved')
'''        
                

# This function returns a fitted covariance matrix. The formula is
# FittedMatrix = LambdaMatrix * PhiMatrix * LambdaMatrix '(transpose) +
# DelMatrix.
# For a more rigorous expression of the formula, see the following paper.
# Jo¨reskog, K. G. (1971). Statistical analysis of sets of congeneric tests.
# Psychometrics, 36 (2), 109-133.
def get_fitcovs(Lambdas, Deltas, PhiMatrix):
    # In LambdaMatrix, the observation variable i has the free parameter lambda
    # if it belongs to the latent variable j, and zero otherwise. Here, all
    # observation variables belong to only one latent variable.
    LambdaMatrix = np.zeros(shape=(NumOfObsVar, NumOfLtnVar))
    for i in range(NumOfObsVar):
        for j in range(NumOfLtnVar):
            if ObsToLtn[i] == j:
                LambdaMatrix[i][j] = Lambdas[i]
    # The one-dimensional array Deltas is converted to the two-dimensional
    # array PhiMatrix, which is a matrix with all zeros except diagonals
    DelMatrix = np.diag(Deltas)
    FittedMatrix = np.empty(shape=(NumOfObsVar, NumOfObsVar))
    FittedMatrix = np.matmul(np.matmul(LambdaMatrix, PhiMatrix),
                             np.transpose(LambdaMatrix)) + DelMatrix
    return FittedMatrix


# The ml (ie, maximum likelihood) function is used as the first argument
# to scipy.optimize.minimize (that is, the objective to be minimized),
# and only a single argument x can be passed. x represents Lambdas,
# Deltas, and PhiVector in order, and the one-dimensional array PhiVector
# is converted to the two-dimensional array PhiMatrix.
def ml(x):
    Lambdas = x[: NumOfObsVar]
    Deltas = x[NumOfObsVar: 2 * NumOfObsVar]
    PhiVector = x[2 * NumOfObsVar:]
    PhiMatrix = get_phimatrix(PhiVector)
    FitCovs = get_fitcovs(Lambdas, Deltas, PhiMatrix)
    return get_mle(FitCovs)


# The null model is a model with variances but no covariances. That is,
# the Lamda matrix is a zero matrix and the non-diagonal part
# of the Phi matrix is all zero. AMOS refers to the null model
# as an independence model.
def ml_null(x):
    FitCovs = np.diag(x)
    return get_mle(FitCovs)


# The ml_sub function is a single-dimensional model for obtaining parameters
# to be used for estimating congeneric reliability.
def ml_sub(x):
    SubLambdas = x[:NOPL]
    SubDeltas = x[NOPL:]
    LambdaMatrix = np.zeros(shape=(NOPL, NOPL))
    for i in range(NOPL):
        for j in range(NOPL):
            if ObsToLtn[i] == j:
                LambdaMatrix[i][j] = SubLambdas[i]
    DelMatrix = np.diag(SubDeltas)
    FitCovs = np.matmul(LambdaMatrix, np.transpose(LambdaMatrix)) + DelMatrix
    return get_mle_sub(FitCovs)


# The get_mle function returns the discrepancy (ie, formula 1) of Jo¨reskog
# (1971). If scipy.optimize.minimize selects the wrong direction once,
# further searching can end up with negative determinants. You can "teach"
# scipy.optimize.minimize to take the right direction by giving arbitrarily
# large numbers to the penalty.
def get_mle(FitCovs):
    PENALTY = 10 ** 10
    if np.linalg.det(FitCovs) > 0:  # determinants of the fitted covariances
        D1 = log(np.linalg.det(FitCovs))  # natural log, not log10
    else:
        D1 = PENALTY
    D2 = np.matrix.trace(np.matmul(ObsCovs, np.linalg.inv(FitCovs)))
    D3 = -log(np.linalg.det(ObsCovs))
    D4 = -NumOfObsVar
    return (D1 + D2 + D3 + D4)


def get_mle_sub(FitCovs):
    PENALTY = 10 ** 10
    if np.linalg.det(FitCovs) > 0:  # determinants of the fitted covariances
        D1 = log(np.linalg.det(FitCovs))  # natural log, not log10
    else:
        D1 = PENALTY
    D2 = np.matrix.trace(np.matmul(ObsCovSubset, np.linalg.inv(FitCovs)))
    D3 = -log(np.linalg.det(ObsCovs))
    D4 = -NOPL
    return (D1 + D2 + D3 + D4)


# The get_phimatrix function converts the one-dimensional array PhiVector
# that constituted the argument of x into a two-dimensional array PhiMatrix.
# One of the following methods is used to determine the scale of
# several manifest variables belonging to the same latent variable.
# 1) Fix one of the lambdas (i.e., path coefficients) to one, or
# 2) Fix the variance of the latent variable to one. Since we are interested in
# the correlations between latent variables, we use the second method to make
# the calculation easier. By doing so, PhiMatrix can be interpreted
# as a correlation matrix between latent variables. The diagonal of PhiMatrix
# is one, and the diagonal above and below are symmetrical to each other.
def get_phimatrix(PhiVector):
    PhiMatrix = np.empty(shape=(NumOfLtnVar, NumOfLtnVar))
    CountAll = 0
    global Cutoff
    if FixPhi == -1:  # -1 means that no phi is fixed
        for i in range(NumOfLtnVar):
            PhiMatrix[i][i] = 1.0
            for j in range(i + 1, NumOfLtnVar):
                PhiMatrix[i][j] = PhiMatrix[j][i] = PhiVector[CountAll]
                CountAll += 1
    else:  # if one of phis must be fixed
        # Only one of all Phi is Fixed and the rest is Free. CountAll
        # increases regardless of whether phi is free or fixed, and
        # CountFree increases only when Phi is Free.
        CountFree = 0
        for i in range(NumOfLtnVar):
            PhiMatrix[i][i] = 1.0
            for j in range(i + 1, NumOfLtnVar):
                if CountAll == FixPhi:
                    if PhiMatrix[i][j] > 0:
                        PhiMatrix[i][j] = Cutoff
                    else:
                        PhiMatrix[i][j] = (-1) * Cutoff
                    PhiMatrix[j][i] = PhiMatrix[i][j]
                    CountAll += 1
                else:
                    PhiMatrix[i][j] = PhiMatrix[j][i] = PhiVector[CountFree]
                    CountAll += 1
                    CountFree += 1
    return PhiMatrix


# For each fit index formula, see the next article.
# Hu, L., & Bentler, P.M. (1999). Cutoff criteria for fit indexes in
# covariance structure analysis: Conventional criteria versus new
# alternatives. Structural Equation Modeling, 6, 1-55.
def get_stats(MinSol, ObjFunVal, NumOfParam, isFitNec=True):
    Chi2 = ObjFunVal * (SampleSize - 1)  # Chi-square
    Df = int(NumOfObsVar * (NumOfObsVar + 1) / 2) - NumOfParam
    Chi2DivDf = Chi2 / Df  # AMOS refers to this as CMIN/DF
    PValue = 1 - chi2.cdf(Chi2, Df)  # p-value
    # Make a dictionary of each statistic.
    Dict = dict()
    Dict["Chi2"] = Chi2
    Dict["Df"] = Df
    Dict["Chi2DivDf"] = Chi2DivDf
    Dict["PValue"] = PValue

    if isFitNec:
        Chi2Null = get_chi2nullmodel()
        DfNull = int(NumOfObsVar * (NumOfObsVar - 1) / 2)
        PValueNull = 1 - chi2.cdf(Chi2Null, DfNull)
        Chi2DivDfNull = Chi2Null / DfNull
        # Normed Fit Index or Bentler-Bonnett Index
        NFI = (Chi2Null - Chi2) / Chi2Null
        # Tucker-Lewis Index or Non-Normed Fit Index
        TLI = (Chi2Null / DfNull - Chi2 / Df) / (Chi2Null / DfNull - 1)
        # Comparative Fit Index, The max function is used to ensure
        # that it has a value between 0 and 1.
        CFI = 1 - max(Chi2 - Df, 0) / max(Chi2 - Df, Chi2Null - DfNull, 0)
        # Incremental Fit Index
        IFI = (Chi2Null - Chi2) / (Chi2Null - Df)
        # Root Mean Squared Error of Approximation
        RMSEA = (max((Chi2 - Df) / (SampleSize - 1), 0) / Df) ** 0.5
        # Standardized root mean squared residual
        Params = x_to_param(MinSol)
        FitCovs = get_fitcovs(Params["Lambdas"], Params["Deltas"],
                              Params["PhiMatrix"])
        sum = 0
        for i in range(NumOfObsVar):
            for j in range(i + 1, NumOfObsVar):
                sum += (ObsCovs[i][j] - FitCovs[i][j]) ** 2 \
                       / (ObsCovs[i][i] * ObsCovs[j][j])
        SRMR = (2 * sum / (NumOfObsVar * (NumOfObsVar + 1))) ** 0.5
        # Akaike Information Criterion
        AIC = Chi2 + NumOfObsVar * (NumOfObsVar + 1) - 2 * Df
        # Bayesian Information Criterion
        BIC = Chi2 + log(SampleSize) \
              * (NumOfObsVar * (NumOfObsVar + 1) / 2 - Df)
        # Sample-Size Adjusted BIC
        SABIC = Chi2 + log((SampleSize + 2) / 24) \
                * (NumOfObsVar * (NumOfObsVar + 1) / 2 - Df)
        # Make a dictionary of each statistic and index.
        Dict["DfNull"] = DfNull
        Dict["Chi2Null"] = Chi2Null
        Dict["DfNull"] = DfNull
        Dict["PValueNull"] = PValueNull
        Dict["Chi2DivDfNull"] = Chi2DivDfNull
        Dict["NFI"] = NFI
        Dict["TLI"] = TLI
        Dict["CFI"] = CFI
        Dict["IFI"] = IFI
        Dict["RMSEA"] = RMSEA
        Dict["SRMR"] = SRMR
        Dict["AIC"] = AIC
        Dict["BIC"] = BIC
        Dict["SABIC"] = SABIC
    return Dict


# obtain standardized lambdas, by dividing its standard deviation
def get_stdlambdas(Lambdas, Deltas):
    # standardized lambdas in vector (one-dimensional array)
    StdLambdas = NumOfObsVar * [0]
    for i in range(NumOfObsVar):
        StdLambdas[i] = Lambdas[i] / sqrt(Lambdas[i] ** 2 + Deltas[i])
    return StdLambdas


# obtain factor pattern matrix, converting standardized lambdas into matrix
def get_patnmatrix(StdLambdas):
    PatnMatrix = np.zeros((NumOfObsVar, NumOfLtnVar))
    for i in range(NumOfObsVar):
        for j in range(NumOfLtnVar):
            if ObsToLtn[i] == j:
                PatnMatrix[i][j] = StdLambdas[i]
    return PatnMatrix


# The null model is a model with variances but no covariances. Thus,
# its parameters are only Deltas (ie, errors), and the number of parameters
# is equal to the number of observed variables. The chi-square value
# of the null model is used to calculate multiple fit indices
# such as NFI, TLI, CFI, and IFI.
def get_chi2nullmodel():
    NumOfParam = NumOfObsVar
    InitialGuess = np.ones(NumOfParam)
    ResultNull = minimize(ml_null, InitialGuess)
    Chi2Null = ResultNull.fun * (SampleSize - 1)
    return Chi2Null


# FixPhi is a global variable that indicates which Phi is fixed. However,
# -1 indicates that all phi are free. Often these variables are passed
# as function arguments, but in scipy.optimize.minimize, you can not pass
# arguments other than x, so we used global variables.
# Greater receives the vector position of the pie
# when the phi is greater than the cutoff. In this case,
# the chi-square test is omitted.
def chi2diff_allpairs(MinSol, ObjFunVal, NumOfParam, Cutoff):
    global FixPhi
    Params = x_to_param(MinSol)
    PhiMatrix = Params["PhiMatrix"]
    Stats = get_stats(MinSol, ObjFunVal, NumOfParam, isFitNec=False)
    # compCutoff indicates the result of comparing each correlation with
    # the cutoff value.
    CompCutoff = NumOfComparison * [""]
    Chi2Fix = NumOfComparison * [0]
    CFIFix = NumOfComparison * [0]
    count = 0
    for i in range(NumOfLtnVar):
        for j in range(i + 1, NumOfLtnVar):
            InitialGuess = np.ones(NumOfParam - 1)
            FixPhi = count
            Result = minimize(ml, InitialGuess)
            StatsFix = get_stats(Result.x, Result.fun, NumOfParam - 1)
            Chi2Fix[count] = StatsFix["Chi2"]
            CFIFix[count] = StatsFix["CFI"]

            # NotSig stores the position of the PhiVector when it is
            # fixed with a cutoff and does not show significant
            # chi-square difference with the original model.
                
            if abs(PhiMatrix[i][j]) > Cutoff:
                # The location of the PhiVector with a correlation
                # greater than the cutoff is stored in the Greater
                CompCutoff[count] = "Greater"
                
            # Then, check the chi-square values while keeping one phi
            # for every other pair. If the difference between the original
            # model and the Chi square is not significant, pass
            # the location to NotSig.
            
            elif StatsFix["Chi2"] < Stats["Chi2"] :
                CompCutoff[count] = "NotSig"
            else:
                CompCutoff[count] = "Sig"
            count += 1
    FixPhi = -1

    Dict = dict()
    Dict["CompCutoff"] = CompCutoff
    Dict["Chi2Fix"] = Chi2Fix
    Dict["CFIFix"] = CFIFix
    return Dict


def get_reliability(Lambdas, Deltas, PhiMatrix):
    # This is mainly a module for obtaining correlated factor reliability
    # and subtest reliability.
    global ObsCovSubset
    global NOPL
    SubRel = NumOfLtnVar * [0]  # Subtest reliability
    ConRel = NumOfLtnVar * [0]  # Congeneric reliability
    TauEqRel = NumOfLtnVar * [0]  # Tau-equivalent reliability
    ParRel = NumOfLtnVar * [0]
    LambdaSum = NumOfLtnVar * [0]
    DeltaSum = NumOfLtnVar * [0]
    NumObsPerLtn = NumOfLtnVar * [0]
    TotalDeltaSum = 0
    TotalVarSum = 0
    FittedMatrix = get_fitcovs(Lambdas, Deltas, PhiMatrix)
    for i in range(NumOfObsVar):
        for j in range(NumOfLtnVar):
            if j == ObsToLtn[i]:
                LambdaSum[j] += Lambdas[i]
                DeltaSum[j] += Deltas[i]
                NumObsPerLtn[j] += 1

    for j in range(NumOfLtnVar):
        TotalDeltaSum += DeltaSum[j]
        SubRel[j] = LambdaSum[j] ** 2 / (LambdaSum[j] ** 2 + DeltaSum[j])

    for i in range(NumOfObsVar):
        for j in range(NumOfObsVar):
            TotalVarSum += FittedMatrix[i][j]
    CFRel = 1 - TotalDeltaSum / TotalVarSum  # Correlated factors reliability


    # This is a module to get tau-equivalent and congeneric reliability.
    # First, obtain the observed covariances for each latent variable.
    for k in range(NumOfLtnVar):
        ObsCovSubset = np.ones(shape=[NumObsPerLtn[k], NumObsPerLtn[k]])
        CountRow = 0
        CountColumn = 0
        i = 0
        j = 0
        while i < NumOfObsVar and CountRow < NumObsPerLtn[k]:
            while j < NumOfObsVar and CountColumn < NumObsPerLtn[k]:
                if ObsToLtn[i] == ObsToLtn[j] == k:
                    ObsCovSubset[CountRow][CountColumn] = ObsCovs[i][j]
                    CountColumn += 1
                j += 1
            i += 1
            j = 0
            if CountColumn > 0:
                CountRow += 1
                CountColumn = 0
        # New ml estimates are obtained to obtain lambdas and deltas
        # in a unidimensional model.
        InitialGuess = np.ones(2 * NumObsPerLtn[k])  # one lambda and one delta
        NOPL = NumObsPerLtn[k]
        # I do not know the exact reason, but I have not made an error here
        # unless I chose "Nelder-Mead". Elsewhere, when I chose "Nelder-Mead"
        # I had the wrong result.
        Result = minimize(ml_sub, InitialGuess, method="Nelder-Mead")
        SubLambdas = Result.x[:NOPL]
        SubDeltas = Result.x[NOPL:]
        LambdaSum = 0
        DeltaSum = 0
        for i in range(NOPL):
            LambdaSum += SubLambdas[i]
            DeltaSum += SubDeltas[i]
        ConRel[k] = LambdaSum ** 2 / (LambdaSum ** 2 + DeltaSum)
        # Now obtain Tau-equivalent reliability.
        TotalSum = 0
        VarSum = 0
        for i in range(NOPL):
            for j in range(NOPL):
                TotalSum += ObsCovSubset[i][j]
                if i == j:
                    VarSum += ObsCovSubset[i][j]
        if NOPL > 1 and TotalSum > 0:
            TauEqRel[k] = (NOPL / (NOPL - 1)) * (1 - VarSum / TotalSum)
        else:
            TauEqRel[k] = 0
        # Now obtain Parallel reliability
        CorSum = 0
        count = 0
        CorAver = 0
        for i in range(NOPL):
            for j in range(i+1,NOPL):
                CorSum += ObsCovSubset[i][j] /\
                sqrt(ObsCovSubset[i][i] * ObsCovSubset[j][j])
                count += 1
        CorAver = CorSum / count
        ParRel[k] = NOPL * CorAver / (1 + (NOPL - 1) * CorAver)

    # And returns the calculated reliability coefficients
    Dict = dict()
    Dict["SubRel"] = SubRel
    Dict["ConRel"] = ConRel
    Dict["TauEqRel"] = TauEqRel
    Dict["ParRel"] = ParRel
    Dict["CFRel"] = CFRel
    return Dict


def get_sscor():
    SMALL_NUMBER = 0.001
    BIG_NUMBER = 100
    SSCor = np.zeros((NumOfLtnVar, NumOfLtnVar))  # Scale score correlation
    SEs = np.zeros((NumOfLtnVar, NumOfLtnVar))  # Standard error of SSCor
    Ts = np.zeros((NumOfLtnVar, NumOfLtnVar))  # T-statistic of SSCor
    Ps = np.zeros((NumOfLtnVar, NumOfLtnVar))  # P-value of SSCor
    CovSum = np.zeros((NumOfLtnVar, NumOfLtnVar))
    VarSum = NumOfLtnVar * [0]

    # The variances and covariances of the observed variables are obtained
    # for each latent variables. This is used to obtain
    # the scale score correlations.
    for i in range(NumOfObsVar):
        for j in range(NumOfObsVar):
            if ObsToLtn[i] == ObsToLtn[j]:
                VarSum[ObsToLtn[i]] += ObsCovs[i][j]
            else:
                CovSum[ObsToLtn[i]][ObsToLtn[j]] += ObsCovs[i][j]

    # Obtaining scale score correlation
    for i in range(NumOfLtnVar):
        for j in range(NumOfLtnVar):
            if i == j:
                SSCor[i][j] = 1
            else:
                SSCor[i][j] = CovSum[i][j] / sqrt(VarSum[i] * VarSum[j])
            SEs[i][j] = sqrt((1 - SSCor[i][j] ** 2) / (SampleSize - 2))
            if SEs[i][j] < SMALL_NUMBER:
                Ts[i][j] = BIG_NUMBER
            else:
                Ts[i][j] = SSCor[i][j] / SEs[i][j]
            # degree of freedom = SampleSize - 2
            # It is a two-tailed test, and it is multiplied by two.
            Ps[i][j] = 2 * (1 - t.cdf(Ts[i][j], SampleSize - 2))

    Dict = dict()
    Dict["SSCor"] = SSCor
    Dict["Ps"] = Ps
    return Dict
    return SSCor

# In scipy.optimize.minimize, use x as the name of all variables
# to optimize. Since its meaning can not be intuitively meaningful,
# this function is responsible for converting to the variable names we use,
# such as Lambda, Delta, and Phi.
def x_to_param(MinSol):
    Lambdas = MinSol[:NumOfObsVar]
    Deltas = MinSol[NumOfObsVar:2 * NumOfObsVar]
    PhiVector = MinSol[2 * NumOfObsVar:]
    PhiMatrix = get_phimatrix(PhiVector)
    Dict = dict()
    Dict["Lambdas"] = Lambdas
    Dict["Deltas"] = Deltas
    Dict["PhiVector"] = PhiVector
    Dict["PhiMatrix"] = PhiMatrix
    return Dict

    
GUI()
