function varargout = ImageProcessingApplication(varargin)
% IMAGEPROCESSINGAPPLICATION MATLAB code for ImageProcessingApplication.fig
%      IMAGEPROCESSINGAPPLICATION, by itself, creates a new IMAGEPROCESSINGAPPLICATION or raises the existing
%      singleton*.
%
%      H = IMAGEPROCESSINGAPPLICATION returns the handle to a new IMAGEPROCESSINGAPPLICATION or the handle to
%      the existing singleton*.
%
%      IMAGEPROCESSINGAPPLICATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEPROCESSINGAPPLICATION.M with the given input arguments.
%
%      IMAGEPROCESSINGAPPLICATION('Property','Value',...) creates a new IMAGEPROCESSINGAPPLICATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageProcessingApplication_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageProcessingApplication_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageProcessingApplication

% Last Modified by GUIDE v2.5 17-Jul-2016 17:22:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageProcessingApplication_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageProcessingApplication_OutputFcn, ...
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

% --- Executes just before ImageProcessingApplication is made visible.
function ImageProcessingApplication_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageProcessingApplication (see VARARGIN)

% Choose default command line output for ImageProcessingApplication
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImageProcessingApplication wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImageProcessingApplication_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Открытие изображения.
function btnOpenImage_Callback(hObject, eventdata, handles)
[fname, pname] = uigetfile('*.*', 'All Files (*.*)');
if fname ~=0
    fullname = strcat(pname, fname);
    image = imread(fullname);
    % imageGray = image(:,:,1);
    % A = im2double(imageGray);
    A = im2double(image);
    figure(1)
    imshow(A)
    handles.A = A;
    guidata(gcbo, handles);
    
    % Condition number
    [~,~,p] = size(A);
    
    for i = 1 : p
        
        switch i
            case 1
                disp('RED')
            case 2
                disp('GREEN')
            case 3
                disp('BLUE')
            otherwise
                disp('UNKNOWN')
        end
        
        s = svd(A(:,:,1));
        disp('Евклидова норма:')
        disp( s(1) )
        disp('Норма Фробениуса:')
        norm_A = sqrt(sum(s.^2));
        disp(norm_A)
        disp('Число обусловленности матрицы исходных данных:')
        Ainv = pinv( A(:,:,1) );
        s = svd(Ainv);
        norm_Ainv = sqrt(sum(s.^2));
        H = norm_A * norm_Ainv;
        disp(H)
    end

end


% --- Импульсный шум.
function btnNoiseImpulse_Callback(hObject, eventdata, handles)
d = str2double(get(handles.edDensity, 'String'));
image = handles.A;
[m,n,p] = size(image);
Anoise = zeros(m, n, p);
for rgb = 1 : p
    Anoise(:,:,rgb) = imnoise(image(:,:,rgb), 'salt & pepper', d); % Импульсный шум
end

% ntsRatio = noiseToSignalRatio( handles.A, Anoise );
% disp('Noise to signal ratio:')
% disp(ntsRatio)
% disp('')

figure(2)
imshow(Anoise)
handles.Anoise = Anoise;
guidata(gcbo, handles);


% --- Гауссов шум.
function btnNoiseGaussian_Callback(hObject, eventdata, handles)
Anoise = imnoise(handles.A, 'gaussian'); % Гауссов шум
figure(2)
imshow(Anoise)
handles.Anoise = Anoise;
guidata(gcbo, handles);


% --- Мультипликативный шум.
function btnNoiseSpeckle_Callback(hObject, eventdata, handles)
Anoise = imnoise(handles.A, 'speckle', 0.04); % Мультипликативный шум
figure(2)
imshow(Anoise)
handles.Anoise = Anoise;
guidata(gcbo, handles);


% --- Усредняющий фильтр.
function btnAverage_Callback(hObject, eventdata, handles)
maskSize = str2num(get(handles.edMaskSize, 'String'));
w = ones(maskSize) / maskSize;
disp(w)
f = handles.Anoise;
gr = imfilter(f, w, 'replicate');
figure, imshow(gr, [ ])


% --- Лапласиан.
function btnLaplacian_Callback(hObject, eventdata, handles)
w = fspecial('laplacian', 0);
f = handles.A;
g = f - imfilter(f, w, 'replicate');
figure, imshow(g)


% --- Медианный фильтр.
function btnMedian_Callback(hObject, eventdata, handles)
maskSize = str2num(get(handles.edMaskSize, 'String'));
f = handles.Anoise;
m = maskSize;
n = m;
%g = ordfilt2(f, median(1 : m*n), ones(m, n));

tic % Начало отсчета времени
g = medfilt2(f, [m, n], 'symmetric');
figure, imshow(g)
T = toc; % Конец отсчета

disp('Filtration time (in seconds):')
disp(T)
disp('')



% --- Построение гистограммы.
function btnHist_Callback(hObject, eventdata, handles)
f = handles.A;
figure, imhist(f)


function edMaskSize_Callback(hObject, eventdata, handles)
% hObject    handle to edMaskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMaskSize as text
%        str2double(get(hObject,'String')) returns contents of edMaskSize as a double


% --- Executes during object creation, after setting all properties.
function edMaskSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edMaskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- SVD.
function btnSVD_Callback(hObject, eventdata, handles)
A = handles.Anoise;
% Filter
level = str2double(get(handles.edLevelSVD, 'String'));
maskSize = str2num(get(handles.edMaskSize, 'String'));

tic % Начало отсчета времени
A1 = SVDImageProcessing(A, maskSize, level);
T = toc; % Конец отсчета

disp('Filtration time (in seconds):')
disp(T)
disp('')

% reverse SVD
recoveryImage = mat2gray(A1); % mat2gray - нормирование
figure, imshow(recoveryImage)


% --- SVD-фильтр импульсных шумов.
function btnSVD_impulse_Callback(hObject, eventdata, handles)
A = handles.Anoise;

tic % Начало отсчета времени
recoveryImage = SVDImageProcessing_impulse(A);
%recoveryImage = SVD_BI_IMPULSE_filter_RGB(A);
T = toc; % Конец отсчета

disp('Filtration time (in seconds):')
disp(T)
disp('')

k = noiseToSignalRatio2(handles.A, handles.Anoise, recoveryImage);

disp('Quality ( k = [noise_sig - orig_sig]/[filtered_sig - orig_sig] ):')
disp(k)
disp('')

figure, imshow(recoveryImage)

% handles.Anoise = recoveryImage;
% guidata(gcbo, handles);


function edLevelSVD_Callback(hObject, eventdata, handles)
% hObject    handle to edLevelSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edLevelSVD as text
%        str2double(get(hObject,'String')) returns contents of edLevelSVD as a double


% --- Executes during object creation, after setting all properties.
function edLevelSVD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edLevelSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edDensity_Callback(hObject, eventdata, handles)
% hObject    handle to edDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edDensity as text
%        str2double(get(hObject,'String')) returns contents of edDensity as a double


% --- Executes during object creation, after setting all properties.
function edDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnCloseAll.
function btnCloseAll_Callback(hObject, eventdata, handles)
close all


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
A = handles.A;
% Filter
level = str2double(get(handles.edLevelSVD, 'String'));
maskSize = str2num(get(handles.edMaskSize, 'String'));

[m, n, p] = size(A);
A1 = zeros(m, n, p);
for rgb = 1 : p
    
    A0 = A(:,:,rgb);
    [U, S, V] = svd(A0);
    
    c = 1;
    
    S1 = c * log(1 + S);

    for i = 1 : m
        A1(:,:,rgb) = A1(:,:,rgb) + U(:,i) * S1(i,i) * V(:,i)';
    end
    
end
B = A + mat2gray(A1);
% reverse SVD
recoveryImage = mat2gray(B); % mat2gray - нормирование
figure, imshow(recoveryImage)