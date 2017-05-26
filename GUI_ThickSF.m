% Copyright (c) 2017.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Pakorn Kanchanawong
% Contact: biekp@nus.edu.sg

function varargout = GUI_ThickSF(varargin)
% GUI_THICKSF MATLAB code for GUI_ThickSF.fig
%      GUI_THICKSF, by itself, creates a new GUI_THICKSF or raises the existing
%      singleton*.
%
%      H = GUI_THICKSF returns the handle to a new GUI_THICKSF or the handle to
%      the existing singleton*.
%
%      GUI_THICKSF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_THICKSF.M with the given input arguments.
%
%      GUI_THICKSF('Property','Value',...) creates a new GUI_THICKSF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ThickSF_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ThickSF_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ThickSF

% Last Modified by GUIDE v2.5 19-Feb-2017 10:10:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ThickSF_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ThickSF_OutputFcn, ...
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


% --- Executes just before GUI_ThickSF is made visible.
function GUI_ThickSF_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ThickSF (see VARARGIN)

% Choose default command line output for GUI_ThickSF
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ThickSF wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ThickSF_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in T_load.
function T_load_Callback(hObject, eventdata, handles)
% hObject    handle to T_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
imgpath=imgetfile;
mkdir data;
save data\imgpath.mat imgpath;
OriginImg = imread(imgpath);
if length(size(OriginImg))==3
    OriginImg = rgb2gray(OriginImg);
end
OriginImg = imadjust(im2uint8(OriginImg));
save data\OriginImg.mat OriginImg;
figure(111);
imshow(OriginImg);axis off;
msgbox('Please Check The Image Loaded');




% --- Executes on button press in T_ROI.
function T_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to T_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
R = str2num(get(handles.T_FilterR,'String'));
NofOrientations_FT = str2num(get(handles.T_FilterAngle,'String'));
save data\R.mat R;
save data\NofOrientations_FT.mat NofOrientations_FT;
load data\OriginImg;
close(figure(111));
I = zeros(size(OriginImg,1)+2*R,  size(OriginImg,2)+2*R);
I(R+1:size(I,1)-R,  R+1:size(I,2)-R) = OriginImg;
figure(111);
imshow(mat2gray(I));
figure(111);ROI_Mask = roipoly;
save data\ROI_Mask.mat ROI_Mask;
msgbox('ROI Selected !');
close(figure(111));



% --- Executes on button press in T_LFTOFT.
function T_LFTOFT_Callback(hObject, eventdata, handles)
% hObject    handle to T_LFTOFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(figure(111));
load data\R.mat R;
load data\NofOrientations_FT.mat NofOrientations_FT;
load data\OriginImg.mat;
load data\ROI_Mask;
[H W] = size(OriginImg);
OriginImg_Margin = uint8(zeros(H+R+R,W+R+R));
OriginImg_Margin(R+1:H+R,R+1:W+R) = OriginImg;

[OFT_Img, LFT_Img, LFT_Orientations] = LFT_OFT_mex(double(OriginImg_Margin),double(R),double(NofOrientations_FT),double(ROI_Mask));
save data\OFT_Img.mat OFT_Img;
save data\LFT_Img.mat LFT_Img;
save data\LFT_Orientations.mat LFT_Orientations;
figure(111);
imshow(mat2gray(OFT_Img));axis off;
msgbox('Transformation Done ! Check the Enhanced Image !');


% --- Executes on button press in T_AutoThresh.
function T_AutoThresh_Callback(hObject, eventdata, handles)
% hObject    handle to T_AutoThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
load data\OFT_Img;
load data\OriginImg;
load data\R;

DefaultFactor = 1;
I = mat2gray(OFT_Img);
t = DefaultFactor*graythresh(I);
if t>=1
    t = graythresh(I);
end
BW = im2bw(I,t);
RawSke = bwmorph(BW,'thin',Inf);
[x y] = find(RawSke==1);
set(handles.T_ThreshValue,'string',t);

close(figure(111));
figure(111);
subplot(1,2,1);
imshow(BW);axis off; title('Segmented Image');
subplot(1,2,2);
imshow(mat2gray(OriginImg));hold on;plot(y-R,x-R,'r.');axis off; title('Original Image with Extracted Skeleton');
scrsz = get(0,'ScreenSize');
set(figure(111),'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])
Allpts = [x y];
AllFragments = RawSke;
save data\Allpts.mat Allpts;
save data\RawSke.mat RawSke;
save data\AllFragments.mat AllFragments;
msgbox('Check Segmented Image (Left) and Extracted Skeleton (Right)');



% --- Executes on button press in T_ThreshBtn.
function T_ThreshBtn_Callback(hObject, eventdata, handles)
% hObject    handle to T_ThreshBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t = str2num(get(handles.T_ThreshValue,'String'));
load data\OFT_Img;
load data\OriginImg;
load data\R;

I = mat2gray(OFT_Img);
BW = im2bw(I,t);
RawSke = bwmorph(BW,'thin',Inf);
[x y] = find(RawSke==1);

close(figure(111));
figure(111);
subplot(1,2,1);
imshow(BW);axis off; title('Segmented Image');
subplot(1,2,2);
imshow(mat2gray(OriginImg));hold on;plot(y-R,x-R,'r.');axis off; title('Original Image with Extracted Skeleton'); % the offset of R/4 is determined by trial and error
scrsz = get(0,'ScreenSize');
set(figure(111),'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
msgbox('Check Segmented Image (Left) and Extracted Skeleton (Right)');
Allpts = [x y];
Allpts = [x y];
AllFragments = RawSke;
save data\Allpts.mat Allpts;
save data\RawSke.mat RawSke;
save data\AllFragments.mat AllFragments;


% --- Executes on button press in T_RemoveJunc.
function T_RemoveJunc_Callback(hObject, eventdata, handles)
% hObject    handle to T_RemoveJunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data\RawSke.mat;
load data\OriginImg.mat;
load data\R.mat;
AllFragments = RawSke;
NOptsMargin = R;
DeleteList = [];

AllFragments(1:NOptsMargin,:) = 0;
AllFragments((end-NOptsMargin):end,:) = 0;
AllFragments(:,1:NOptsMargin) = 0;
AllFragments(:,(end-NOptsMargin):end) = 0;

[x y] = find(AllFragments==1);
Allpts = [x y];
% remove margin information above

% remove crossing points, a N-by-N region located at the base point will be removed
R_Junc = floor((str2num(get(handles.T_JuncSize,'String'))-1)/2);
Size_Junc = str2num(get(handles.T_JuncSize,'String'));
N1 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)-1));
N2 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)));
N3 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)+1));
N4 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1),Allpts(:,2)-1));
N5 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1),Allpts(:,2)+1));
N6 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)-1));
N7 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)));
N8 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)+1));

N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;

CrPts = Allpts(find(N>=3),1:2);

h = waitbar(0,'Removing crossing points');
for i = 1:size(CrPts,1)
    waitbar(i/size(CrPts,1),h);
    AllFragments(CrPts(i,1)-R_Junc:CrPts(i,1)+R_Junc,CrPts(i,2)-R_Junc:CrPts(i,2)+R_Junc) = 0;
end

close(h);
[x y] = find(AllFragments==1);
Allpts = [x y];

N1 = AllFragments(sub2ind(size(AllFragments),x-1,y-1));
N2 = AllFragments(sub2ind(size(AllFragments),x-1,y));
N3 = AllFragments(sub2ind(size(AllFragments),x-1,y+1));
N4 = AllFragments(sub2ind(size(AllFragments),x,y-1));
N5 = AllFragments(sub2ind(size(AllFragments),x,y+1));
N6 = AllFragments(sub2ind(size(AllFragments),x+1,y-1));
N7 = AllFragments(sub2ind(size(AllFragments),x+1,y));
N8 = AllFragments(sub2ind(size(AllFragments),x+1,y+1));

N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;

SinglePts = Allpts(find(N==0),1:2);
% remove single points
h = waitbar(0,'Removing single points');
for i = 1:size(SinglePts,1)
    waitbar(i/size(SinglePts,1),h);
    AllFragments(SinglePts(i,1),SinglePts(i,2)) = 0;
end
close(h);
[x y] = find(AllFragments==1);
Allpts = [x y];
[L num] = bwlabel(AllFragments,8);

RawCrPts = CrPts;
save data\L.mat L num;
save data\AllFragments.mat AllFragments;
save data\Allpts.mat Allpts;
save data\Size_Junc.mat Size_Junc;
save data\RawCrPts.mat RawCrPts;
close(figure(111));
figure(111);
imshow(mat2gray(OriginImg));hold on;plot(y-R,x-R,'r.');axis off;
% remove clusters of single points


% --- Executes on button press in T_ShortFragRemove.
function T_ShortFragRemove_Callback(hObject, eventdata, handles)
% hObject    handle to T_ShortFragRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data\AllFragments.mat;
load data\R.mat
load data\OriginImg.mat
load data\Allpts.mat;
load data\L.mat;
MIN_FragmentLength = str2num(get(handles.T_ShortFrag,'String'));
AllFragments = bwareaopen(AllFragments, MIN_FragmentLength);
[L num] = bwlabel(AllFragments,8);
[x y] = find(AllFragments==1);
Allpts = [x y];
save data\AllFragments.mat AllFragments;
save data\Allpts.mat Allpts;
save data\L.mat L num;
close(figure(111));
figure(111);
imshow(mat2gray(OriginImg));hold on;
plot(y-R,x-R,'r.');axis off;
% remove very short fragment


% --- Executes on button press in T_IterBtn.
function T_IterBtn_Callback(hObject, eventdata, handles)
% hObject    handle to T_IterBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data\OriginImg.mat;
load data\AllFragments.mat;
load data\OFT_Img.mat;
load data\R.mat;
load data\ROI_Mask.mat;
load data\NofOrientations_FT.mat;

Thresh = str2num(get(handles.T_ThreshValue,'String'));
JuncSize = str2num(get(handles.T_JuncSize,'String'));
MIN_FragmentLength = str2num(get(handles.T_ShortFrag,'String'));
iteration = get(handles.T_IterN,'Value');
Iter_RemoveR = 3; % unit: pixels
R_Junc = (JuncSize-1)/2;
AllFragments = IterGenFragment(OriginImg, ...
    AllFragments, ...
    iteration, ...
    R, ...
    ROI_Mask, ...
    NofOrientations_FT, ...
    Iter_RemoveR, ...
    Thresh, ...
    R_Junc, ...
    MIN_FragmentLength);
[L num] = bwlabel(AllFragments);
[x y] = find(AllFragments==1);
Allpts = [x y];
close(111);
figure(111);
imshow(mat2gray(OriginImg));hold on;
plot(y-R,x-R,'r.');axis off;

save data\AllFragments.mat AllFragments;
save data\Allpts.mat Allpts;
save data\L.mat L num;


% --- Executes on button press in T_TipReg.
function T_TipReg_Callback(hObject, eventdata, handles)
% hObject    handle to T_TipReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data\AllFragments.mat;
load data\Allpts.mat;
load data\L.mat;
load data\R.mat;
LL = L;
SkeTermi = bwmorph(AllFragments,'endpoints');
[x y] = find(SkeTermi==1);
all_tips = [x y];

all_tips(:,3) = L(sub2ind(size(L),all_tips(:,1),all_tips(:,2)));

% may register tip orientation using parallel computing below
RegR = 2*R; % to register the direction, we need to consider only a local region around a tip
save data\RegR.mat RegR;
tempInfo = zeros(size(all_tips,1),3);
h = waitbar(0,'Registering Tips ...');
for i = 1:size(all_tips,1)
    waitbar(i/size(all_tips,1),h);
    tempInfo(i,:) = TipReg(LL, all_tips(i,3), all_tips(i,1:2), RegR);
end
close(h);

all_tips = [all_tips, tempInfo];
% may register tip orientation using parallel computing above

close(figure(111));
figure(111);
imshow(mat2gray(AllFragments));hold on; axis off;plot(all_tips(:,2),all_tips(:,1),'r+');
for i = 1:100
    idx = ceil(rand * size(all_tips,1));
    text(all_tips(idx,2),all_tips(idx,1),[num2str(all_tips(idx,4))],'color','g');
end
save data\all_tips.mat all_tips;


% --- Executes on button press in T_SearchPrev.
function T_SearchPrev_Callback(hObject, eventdata, handles)
% hObject    handle to T_SearchPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(111));
load data\AllFragments;
load data\all_tips;
figure(111);
imshow(mat2gray(AllFragments));axis off;title('Please double-click the tip where you want to check the region for searching');
[Y1 X1]= ginput(1);
k = dsearchn(all_tips(:,1:2),[X1 Y1]);
k = k(1);
BasePt = all_tips(k,:);
FanR = str2num(get(handles.T_SearchR,'String'));
FanAngle = str2num(get(handles.T_SearchAngle,'String'));
FanEdge = FanPreview(FanAngle, AllFragments, FanR, BasePt(4), BasePt(1:2));
close(figure(111));
figure(111);
imshow(mat2gray(AllFragments));hold on; axis off;
plot(FanEdge(:,2), FanEdge(:,1), 'g', 'LineWidth', 2);



% --- Executes on button press in T_TipPairing.
function T_TipPairing_Callback(hObject, eventdata, handles)
% hObject    handle to T_TipPairing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(figure(111));
load data\all_tips.mat;
load data\AllFragments;
load data\L.mat;

FanR = str2num(get(handles.T_SearchR,'String'));
FanAngle = str2num(get(handles.T_SearchAngle,'String'));

MIN_Angle_Diff = str2num(get(handles.T_DirDiff,'String'));
MIN_Dist = str2num(get(handles.T_SearchR,'String'));
set(handles.T_ShortFila,'string',MIN_Dist);
MIN_GapAngle_Diff = str2num(get(handles.T_GapDiff,'String'));

MIN_Info = [MIN_Angle_Diff,MIN_Dist,MIN_GapAngle_Diff];

label_list = all_tips(:,3); % will be updated; help to find another tip of non-first filament during searching
L_GlobalIndex = zeros(size(L));
all_tips(:,7) = (1:size(all_tips,1))';
L_GlobalIndex(sub2ind(size(L_GlobalIndex),all_tips(:,1),all_tips(:,2))) = all_tips(:,7);
% go through all tips and register all searched information

new_partner_list = zeros(size(all_tips,1),100); % reserve the memory; a maximum of 100 partner tips allowed
C1weight = str2num(get(handles.T_DirDiffW,'String'));
C3weight = str2num(get(handles.T_GapDiffW,'String'));

h = waitbar(0,'Tip Searching in Progress ...');
for i = 1:size(all_tips,1)
    waitbar(i/size(all_tips,1),h);
    new_partner_list(i,:) = local_search(L_GlobalIndex, ...            % image of skeleton and its tips are labeled with global index
        label_list, ...               % this list will be frequently updated (if a tip has been used, it will be removed from this list)
        FanAngle, ...                 % size of the fan-shape searching region
        all_tips, ...                 % information of all tips
        FanR,  ...                    % this radius defines the range for searching
        all_tips(i,4), ...            % orientation of that current tip
        all_tips(i,1:2), ...          % coordinate of current tip
        MIN_Info, ...                 % Minimum conditions that should satisfy
        C1weight, ...                 % weight of croterion 1
        C3weight);                    % weight of croterion 3
    
end
close(h);

for i = 1:size(new_partner_list,2)
    if sum(new_partner_list(:,i))==0
        new_partner_list(:,i:end) = [];
        break;
    end
end
save data\new_partner_list.mat new_partner_list;

all_tips(:,7) = 1:size(all_tips,1); % add global index to all tips
all_tips(:,8) = ones(size(all_tips,1),1); % this reserves to indicates number of lives

h = waitbar(0,'Add Number of Lives for Each Filament ...');
% add the number of lives & partner tips to all tips
all_tips(:,9:end) = []; % clear previous record
for i = 1:size(all_tips,1)
    waitbar(i/size(all_tips,1),h);
    temp_partner = new_partner_list(i,:);
    temp_partner(find(temp_partner==0)) = [];
    if ~isempty(temp_partner)
        all_tips(i,8) = length(temp_partner);
        all_tips(i,9:(8 + length(temp_partner))) = temp_partner;
    end
end
close(h);
% add the number of lives & partner tips to all tips
[L num] = bwlabel(AllFragments,8);
L_GlobalIndex = zeros(size(L));
L_GlobalIndex(sub2ind(size(L),all_tips(:,1),all_tips(:,2))) = 1:length(all_tips(:,1));
save data\L_GlobalIndex.mat L_GlobalIndex;
msgbox('Tip Search Done ! Please Proceed to GROUPING !');

% assign number of lives to each fragments according to the max number of lives of its two tips
for i = 1:num
    pair_index = find(all_tips(:,3)==i);
    if isempty(pair_index)
        continue;
    end
    all_tips(pair_index(1),8) = max(all_tips(pair_index,8));
    all_tips(pair_index(2),8) = max(all_tips(pair_index,8));
end

all_tips = biDirPairing(all_tips);

save data\all_tips.mat all_tips;

% ****** structure of all_tips so far ******   #: Number
%           colume1            colume2           colume3         colume4         colume5            colume6        colume7          colume8          colume9               colume10                colume11          ...
% tip1     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% tip2     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% tip3     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% tip4     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% ......
% ****** structure of all_tips so far ******   #: Number
% assign number of lives to each fragments according to the max number of lives of its two tips



% --- Executes on button press in T_Grp.
function T_Grp_Callback(hObject, eventdata, handles)
% hObject    handle to T_Grp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% filamentous fragment grouping starts
close(figure(111));
load data\L.mat;
load data\all_tips.mat;

NofLives = all_tips(:,8);
NofLives(find(NofLives==0)) = 1;
all_tips(:,8) = NofLives;
% assign number of lives to each fragments according to the max number of lives of its two tips

new_partner_list = all_tips(:,9:end); % update the list of partners

all_filament = [];
all_connects = [];
[x y] = find(L~=0);
Lpts = [x, y, L(sub2ind(size(L),x,y))];
h = waitbar(0,'Grouping in Progress');
tic;
for i = 1:num
    waitbar(i/num,h);
    L_list = [];
    label_list = all_tips(:,3); % labels (corresponding to L) of all fragments the global index of which corresponds to all_tips
    xx = Lpts(find(Lpts(:,3)==i),1);
    yy = Lpts(find(Lpts(:,3)==i),2);
    if ~isempty(xx) % if this fragment still exists
        tips_index = find(all_tips(:,3)==i);% find global index
        newSinglefilament = []; % used to store a single filament
        newSinglefilament = [newSinglefilament; [xx yy]]; % add up to the single filament
        L_list = [L_list  i];
        tip_dir1 = tips_index(1); % change to the first tip
        tip_dir2 = tips_index(2);
        connect_tips = all_tips(tip_dir1,1:2); % the design of this connect_tips is for quick filling between gaps later
        if NofLives(tip_dir1)<2
            Lpts(find(Lpts(:,3)==i),3)=0;
            label_list(tips_index(1)) = 0;
            label_list(tips_index(2)) = 0;
            new_partner_list(find(new_partner_list==tip_dir1)) = 0;
            new_partner_list(find(new_partner_list==tip_dir2)) = 0;
        else
            NofLives(tips_index(1)) = NofLives(tips_index(1)) - 1; % reduce its number of lives by 1 if it is used
            NofLives(tips_index(2)) = NofLives(tips_index(2)) - 1; % reduce its number of lives by 1 if it is used
        end
        while 1+1==2
            can_index = new_partner_list(tip_dir1,:); % get global indices of possible partner tips
            can_index(find(can_index==0)) = []; % remove zeros
            if isempty(can_index) || ~isempty(find(L_list==all_tips(can_index(1),3)))  % if the current tip has no possible partner tip
                break;
            else
                optimal_index = can_index(1);
                xx = Lpts(find(Lpts(:,3)==label_list(optimal_index)),1);
                yy = Lpts(find(Lpts(:,3)==label_list(optimal_index)),2);
                newSinglefilament = [newSinglefilament; [xx yy]]; % add new fragment to the current single filament
                L_list = [L_list  label_list(optimal_index)];
                connect_tips = [connect_tips;all_tips(optimal_index,1:2)];
                ll = NofLives(optimal_index);
                if NofLives(optimal_index)<2
                    Lpts(find(Lpts(:,3)==label_list(optimal_index)),3)=0;
                    new_partner_list(find(new_partner_list==tip_dir1)) = 0;
                    label_list(optimal_index) = 0; % remove this label but leave the other one
                    new_partner_list(find(new_partner_list==optimal_index)) = 0;
                else
                    NofLives(find(label_list==all_tips(optimal_index,3))) = NofLives(find(label_list==all_tips(optimal_index,3))) - 1; % reduce its number of lives by 1 if it is used
                end
                new_partner_list(tip_dir1,  find(new_partner_list(tip_dir1,:)==optimal_index)  ) = 0; % remove the tip from the current list of partner tips
                current_label = all_tips(optimal_index,3); % get the current label
                tip_dir1 = find(label_list==current_label); % find the other tip since only one left now
                if length(tip_dir1)==2
                    tip_dir1 = tip_dir1(find(tip_dir1~=optimal_index)); % find the other end of the fragment
                end
                if ll<2
                    new_partner_list(find(new_partner_list==tip_dir1)) = 0;
                    label_list(tip_dir1) = 0;
                end
                connect_tips = [connect_tips;all_tips(tip_dir1,1:2)];
            end
        end
        tip_dir2 = tips_index(2);
        connect_tips = [all_tips(tip_dir2,1:2);connect_tips];
        while 1+1==2
            can_index = new_partner_list(tip_dir2,:); % get global indices of possible partner tips
            can_index(find(can_index==0)) = []; % remove zeros
            if isempty(can_index) || ~isempty(find(L_list==all_tips(can_index(1),3)))% if the current tip has no possible partner tip
                break;
            else
                optimal_index = can_index(1);
                xx = Lpts(find(Lpts(:,3)==label_list(optimal_index)),1);
                yy = Lpts(find(Lpts(:,3)==label_list(optimal_index)),2);
                newSinglefilament = [newSinglefilament; [xx yy]]; % add new fragment to the current single filament
                L_list = [L_list  label_list(optimal_index)];
                connect_tips = [all_tips(optimal_index,1:2);connect_tips];
                ll = NofLives(optimal_index);
                if NofLives(optimal_index)<2
                    Lpts(find(Lpts(:,3)==label_list(optimal_index)),3)=0;
                    new_partner_list(find(new_partner_list==tip_dir2)) = 0;
                    label_list(optimal_index) = 0; % remove this label but leave the other one
                    new_partner_list(find(new_partner_list==optimal_index)) = 0;
                else
                    NofLives(find(label_list==all_tips(optimal_index,3))) = NofLives(find(label_list==all_tips(optimal_index,3))) - 1; % reduce its number of lives by 1 if it is used
                end
                new_partner_list(tip_dir2,  find(new_partner_list(tip_dir2,:)==optimal_index)  ) = 0; % remove the tip from the current list of partner tips
                current_label = all_tips(optimal_index,3); % get the current label
                tip_dir2 = find(label_list==current_label); % find the other tip since only one left now
                if length(tip_dir2)==2
                    tip_dir2 = tip_dir2(find(tip_dir2~=optimal_index)); % find to other end of the fragment
                end
                if ll<2
                    new_partner_list(find(new_partner_list==tip_dir2)) = 0;
                    label_list(tip_dir2) = 0;
                end
                connect_tips = [all_tips(tip_dir2,1:2);connect_tips];
            end
        end
        if isempty(all_filament) && isempty(all_connects)
            all_filament(1:size(newSinglefilament,1),1:2,1) = newSinglefilament;
            all_connects(1:size(connect_tips,1),1:2,1) = connect_tips;
        else
            all_filament(1:size(newSinglefilament,1),1:2,size(all_filament,3)+1) = newSinglefilament;
            all_connects(1:size(connect_tips,1),1:2,size(all_connects,3)+1) = connect_tips;
        end
    end
end
toc;
msgbox('Grouping Done !');
close(h);
close(figure(111));
% filamentous fragment grouping ends
save data\all_filament.mat all_filament;
save data\all_connects.mat all_connects;
all_connects_shortlist = all_connects;
save data\all_connects_shortlist.mat all_connects_shortlist;



% --- Executes on button press in T_Sort.
function T_Sort_Callback(hObject, eventdata, handles)
% hObject    handle to T_Sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(figure(111));
load data\all_filament.mat;
load data\all_connects.mat;
load data\L;

maskI = zeros(size(L));
Fullength = 5*size(all_filament,1);
all_sorted_filament = zeros(Fullength,2,size(all_filament,3));

h = waitbar(0,'Analysis in Progress ...');
for i = 1:size(all_filament,3)
    waitbar(i/size(all_filament,3),h);
    [all_sorted_filament(:,:,i)  AnalysisInfo(i,:)]= SortFilament(all_filament(:,:,i), maskI, all_connects(:,:,i), Fullength);
end
close(h);

for i = 1:size(all_sorted_filament,1)
    temp = all_sorted_filament(i,:,:);
    if(sum(temp(:)))==0
        all_sorted_filament(i:end,:,:) = [];
        break;
    end
end
msgbox('Analysis and Sorting Done !');

save data\all_sorted_filament.mat all_sorted_filament;
save data\AnalysisInfo.mat AnalysisInfo;



% --- Executes on button press in T_ShortFilaRemove.
function T_ShortFilaRemove_Callback(hObject, eventdata, handles)
% hObject    handle to T_ShortFilaRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(figure(111));

load data\AnalysisInfo;
load data\all_sorted_filament;
load data\AllFragments.mat;
load data\all_connects.mat;
load data\OriginImg;
load data\R;

ShortFilament= str2num(get(handles.T_ShortFila,'String'));

RemoveIdx = [];

h = waitbar(0,'Removing Short Filaments...');
for i = 1:size(all_sorted_filament,3)
    waitbar(i/size(all_sorted_filament,3),h);
    if AnalysisInfo(i,2) <= ShortFilament
        RemoveIdx = [RemoveIdx i];
    end
end
close(h);
all_sorted_filament(:,:,RemoveIdx) = [];
AnalysisInfo(RemoveIdx,:) = [];
all_connects(:,:,RemoveIdx) = [];

% remove ungrouped below
RemoveUngrp = get(handles.T_ungrp,'Value');
if RemoveUngrp==1
    AllFragments = im2bw(AllFragments);
    RemoveIdx = [];
    h = waitbar(0,'Removing Ungrouped Fragments...');
    for i = 1:size(all_sorted_filament,3)
        waitbar(i/size(all_sorted_filament,3),h);
        temp = all_sorted_filament(:,:,i);
        temp(find(temp(:,1)==0),:) = [];
        if size(temp,1)==sum(AllFragments(sub2ind(size(AllFragments),temp(:,1),temp(:,2))))
            RemoveIdx = [RemoveIdx  i];
        end
    end
    close(h);
    all_sorted_filament(:,:,RemoveIdx) = [];
    AnalysisInfo(RemoveIdx,:) = [];
    all_connects(:,:,RemoveIdx) = [];
end
% remove ungrouped above
all_connects_shortlist = all_connects;
save data\AnalysisInfo.mat AnalysisInfo;
save data\all_sorted_filament.mat all_sorted_filament;
save data\all_connects_shortlist.mat all_connects_shortlist;

figure(111);
imshow(OriginImg);hold on;
for i = 1:size(all_sorted_filament,3)
    x = all_sorted_filament(:,1,i);x(x==0) = [];
    y = all_sorted_filament(:,2,i);y(y==0) = [];
    figure(111);plot(y-R,x-R,'color',rand(1,3),'LineWidth',2);
end
    




% --- Executes on button press in T_MeasureExport.
function T_MeasureExport_Callback(hObject, eventdata, handles)
% hObject    handle to T_MeasureExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;

load data\OFT_Img;
load data\R;
load data\OriginImg;
load data\all_sorted_filament;
load data\AnalysisInfo;
close(figure(111));
[H W] = size(OriginImg); 
new = uint8(zeros(H+R+R,W+R+R));
new((R+1):(R+H),(R+1):(R+W)) = OriginImg;
OriginImg = new;

[Gmag,Gdir] = imgradient(OriginImg,'sobel');

AllS_lower = str2num(get(handles.T_WidthLow,'String'));
AllS_upper = str2num(get(handles.T_WidthUp,'String'));
figure(111);imshow(OriginImg);hold on;axis off;
for i = 1:size(all_sorted_filament,3)
    temp = all_sorted_filament(:,:,i);
    temp(find(temp(:,1)==0),:) = [];
    D = zeros(size(OFT_Img));
    D(sub2ind(size(D),temp(:,1),temp(:,2))) = 1;
    D = bwdist(D);
    AllS = AllS_lower:0.2:AllS_upper;
    AllG = [];
    for j = 1:length(AllS)
        [x y] = find(D<AllS(j));
        SingleF = zeros(size(D));
        SingleF(sub2ind(size(SingleF),x,y)) = 1;
        B=bwboundaries(SingleF);
        B = B{1};
        G = mean(Gmag(sub2ind(size(Gmag),B(:,1),B(:,2))));
        AllG(j) = G;
    end
    idx = find(AllG==max(AllG));
    idx = idx(1);
    S = AllS(idx);
    [x y] = find(D<AllS(idx));
    SingleF = zeros(size(D));
    SingleF(sub2ind(size(SingleF),x,y)) = 1;
    B=bwboundaries(SingleF);
    B = B{1};
    figure(111);plot(B(:,2),B(:,1),'color','r');
    display(['Processing Filament ',num2str(i),'/',num2str(size(all_sorted_filament,3)),'. Size = ',num2str(2*S),' pixels']);
    FinalS(i) = 2*S;
end

for i = 1:size(all_sorted_filament,3)
    temp = all_sorted_filament(:,:,i);
    temp(find(temp(:,1)==0),:) = [];
    figure(111);
    text(temp(1,2)-R,temp(1,1)-R,['(',num2str(i),')'],'color','g');
end

mkdir result;

titles = {'Filament ID','Orientation','Total Length','End-to-End Distance','Centroid X','Centroid Y','Width'};
titles = [titles,repmat({'X','Y'},[1 size(all_sorted_filament,1)]) ];
InfoExcel = zeros(size(all_sorted_filament,3), size(all_sorted_filament,1)*2+6);

for i = 1:size(all_sorted_filament,3)
    curr_filament = all_sorted_filament(:,:,i);
    curr_filament(find(curr_filament(:,1)==0),:) = [];
    InfoExcel(i, 1:6) = [AnalysisInfo(i,:)  FinalS(i)];
    temp = zeros(1,size(all_sorted_filament,1)*2);
    temp((1:size(all_sorted_filament,1))*2-1) = all_sorted_filament(:,1,i);
    temp((1:size(all_sorted_filament,1))*2) = all_sorted_filament(:,2,i);
    InfoExcel(i, 7:end) = temp;
end
FragID = 1:size(all_sorted_filament,3);
FragID = FragID(:);
ThickSF_Info = [FragID, InfoExcel];
save data\ThickSF_Info.mat ThickSF_Info;
InfoExcel = [titles;num2cell(ThickSF_Info)];

xlswrite('result\ThickSF_Info.xlsx',InfoExcel,1,'A1');
msgbox('Stress Fiber Measurement and Data Exportation Done !');





function T_FilterR_Callback(hObject, eventdata, handles)
% hObject    handle to T_FilterR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_FilterR as text
%        str2double(get(hObject,'String')) returns contents of T_FilterR as a double


% --- Executes during object creation, after setting all properties.
function T_FilterR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_FilterR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_FilterAngle_Callback(hObject, eventdata, handles)
% hObject    handle to T_FilterAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_FilterAngle as text
%        str2double(get(hObject,'String')) returns contents of T_FilterAngle as a double


% --- Executes during object creation, after setting all properties.
function T_FilterAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_FilterAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_ThreshValue_Callback(hObject, eventdata, handles)
% hObject    handle to T_ThreshValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_ThreshValue as text
%        str2double(get(hObject,'String')) returns contents of T_ThreshValue as a double


% --- Executes during object creation, after setting all properties.
function T_ThreshValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_ThreshValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_JuncSize_Callback(hObject, eventdata, handles)
% hObject    handle to T_JuncSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_JuncSize as text
%        str2double(get(hObject,'String')) returns contents of T_JuncSize as a double


% --- Executes during object creation, after setting all properties.
function T_JuncSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_JuncSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_ShortFrag_Callback(hObject, eventdata, handles)
% hObject    handle to T_ShortFrag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_ShortFrag as text
%        str2double(get(hObject,'String')) returns contents of T_ShortFrag as a double


% --- Executes during object creation, after setting all properties.
function T_ShortFrag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_ShortFrag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_SearchR_Callback(hObject, eventdata, handles)
% hObject    handle to T_SearchR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_SearchR as text
%        str2double(get(hObject,'String')) returns contents of T_SearchR as a double


% --- Executes during object creation, after setting all properties.
function T_SearchR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_SearchR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_SearchAngle_Callback(hObject, eventdata, handles)
% hObject    handle to T_SearchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_SearchAngle as text
%        str2double(get(hObject,'String')) returns contents of T_SearchAngle as a double


% --- Executes during object creation, after setting all properties.
function T_SearchAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_SearchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_DirDiff_Callback(hObject, eventdata, handles)
% hObject    handle to T_DirDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_DirDiff as text
%        str2double(get(hObject,'String')) returns contents of T_DirDiff as a double


% --- Executes during object creation, after setting all properties.
function T_DirDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_DirDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_GapDiff_Callback(hObject, eventdata, handles)
% hObject    handle to T_GapDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_GapDiff as text
%        str2double(get(hObject,'String')) returns contents of T_GapDiff as a double


% --- Executes during object creation, after setting all properties.
function T_GapDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_GapDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_DirDiffW_Callback(hObject, eventdata, handles)
% hObject    handle to T_DirDiffW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_DirDiffW as text
%        str2double(get(hObject,'String')) returns contents of T_DirDiffW as a double


% --- Executes during object creation, after setting all properties.
function T_DirDiffW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_DirDiffW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_GapDiffW_Callback(hObject, eventdata, handles)
% hObject    handle to T_GapDiffW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_GapDiffW as text
%        str2double(get(hObject,'String')) returns contents of T_GapDiffW as a double


% --- Executes during object creation, after setting all properties.
function T_GapDiffW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_GapDiffW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_ShortFila_Callback(hObject, eventdata, handles)
% hObject    handle to T_ShortFila (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_ShortFila as text
%        str2double(get(hObject,'String')) returns contents of T_ShortFila as a double


% --- Executes during object creation, after setting all properties.
function T_ShortFila_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_ShortFila (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_WidthLow_Callback(hObject, eventdata, handles)
% hObject    handle to T_WidthLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_WidthLow as text
%        str2double(get(hObject,'String')) returns contents of T_WidthLow as a double


% --- Executes during object creation, after setting all properties.
function T_WidthLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_WidthLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in T_IterN.
function T_IterN_Callback(hObject, eventdata, handles)
% hObject    handle to T_IterN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns T_IterN contents as cell array
%        contents{get(hObject,'Value')} returns selected item from T_IterN


% --- Executes during object creation, after setting all properties.
function T_IterN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_IterN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in T_ungrp.
function T_ungrp_Callback(hObject, eventdata, handles)
% hObject    handle to T_ungrp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T_ungrp



function T_WidthUp_Callback(hObject, eventdata, handles)
% hObject    handle to T_WidthUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_WidthUp as text
%        str2double(get(hObject,'String')) returns contents of T_WidthUp as a double


% --- Executes during object creation, after setting all properties.
function T_WidthUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_WidthUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in T_BackBtn.
function T_BackBtn_Callback(hObject, eventdata, handles)
% hObject    handle to T_BackBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;close all;
SFEX;
