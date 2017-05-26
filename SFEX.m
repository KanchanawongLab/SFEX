% Copyright (c) 2017.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Pakorn Kanchanawong
% Contact: biekp@nus.edu.sg

function varargout = SFEX(varargin)
% SFEX MATLAB code for SFEX.fig
%      SFEX, by itself, creates a new SFEX or raises the existing
%      singleton*.
%
%      H = SFEX returns the handle to a new SFEX or the handle to
%      the existing singleton*.
%
%      SFEX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SFEX.M with the given input arguments.
%
%      SFEX('Property','Value',...) creates a new SFEX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SFEX_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SFEX_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SFEX

% Last Modified by GUIDE v2.5 18-Feb-2017 11:58:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SFEX_OpeningFcn, ...
                   'gui_OutputFcn',  @SFEX_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SFEX is made visible.
function SFEX_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SFEX (see VARARGIN)

% Choose default command line output for SFEX
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SFEX wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SFEX_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in toGUI_ThickSF.
function toGUI_ThickSF_Callback(hObject, eventdata, handles)
% hObject    handle to toGUI_ThickSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;close all;
GUI_ThickSF;


% --- Executes on button press in toGUI_SFnetwork.
function toGUI_SFnetwork_Callback(hObject, eventdata, handles)
% hObject    handle to toGUI_SFnetwork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;close all;
GUI_SFnetwork;
