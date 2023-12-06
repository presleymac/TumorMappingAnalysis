close all 
clear all 
clc 

matlab_folder = pwd;%remembers the current folder so it can come back to it
%load required variables

data_folder = uigetdir('','Select folder with image files');   

cd(data_folder); 
answer = inputdlg('How many hours from now should files be converted?');
disp('delayed for')
disp(strcat(cell2mat(answer), ' hours'))

pause(str2num(cell2mat(answer))*60*60);

files = dir('*.czi');    

%%
for int_file = 1:size(files,1)   
    
    [~,shortfile] = fileparts(files(int_file).name);    
    mkdir(shortfile);                            
    cd(shortfile);                               
    
    origfilename_ch1 = strcat(shortfile, '_iso_ch1.tif'); 
    origfilename_ch2 = strcat(shortfile, '_iso_ch2.tif');
    origfilename_ch3 = strcat(shortfile, '_iso_ch3.tif');
    origfilename_ch4 = strcat(shortfile, '_iso_ch4.tif');
    origfilename_ch5 = strcat(shortfile, '_iso_ch5.tif');
    filenames = [origfilename_ch1; origfilename_ch2; origfilename_ch3; origfilename_ch4; origfilename_ch5];  
    
    metadatafilename = strcat(shortfile, '_iso_info.csv');
    metadatafullfilename = strcat(shortfile, '_iso_fullmeta.csv'); 
    if exist(metadatafullfilename,'file')==2    
        disp('skipping');
        disp(origfilename_ch1);
    else
        
        
        tic
        cd(data_folder);                   
        data = bfopen(files(int_file).name);       
        disp('reading time');
        toc
        cd(shortfile);
        omeMeta = data{1, 4};
        stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
        stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
        stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
        numchannels = omeMeta.getPixelsSizeC(0).getValue();
        timepoints = omeMeta.getPixelsSizeT(0).getValue();
        
        metadata = data{1, 2};
        metadataKeys = metadata.keySet().iterator();
        metadataarray = {};
        for m=1:metadata.size()
            key = metadataKeys.nextElement();
            value = metadata.get(key);
            %   fprintf('%s = %s\n', key, value)
            metadataarray = [metadataarray; {key} {value}];
        end
       
        [metadatakeys_cell, I]= sort( metadataarray(:,1));
        metadatavalues_cell =  metadataarray(I,2);
        
        
        T2 = table(metadatakeys_cell,metadatavalues_cell);
        
        [laser_wavelength, laser_powers, cam_int_time, filterset, ...
            trackconfig, leftrightalignment, LSthickness, imagezoom, ... 
            objective] = extractMetaData(metadatakeys_cell,metadatavalues_cell);
        %         strphysSizeX = (metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'));
        %         strphysSizeY = (metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'));
        %         strphysSizeZ = (metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1'));
        
        
        physSizeX = str2num(metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'));
        physSizeY = str2num(metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'));
        physSizeZ = str2num(metadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1'));
        
        newSizeX = round(stackSizeX*physSizeX/physSizeZ);
        newSizeY = round(stackSizeY*physSizeY/physSizeZ);
        newSizeZ = stackSizeZ;
        
        newphysX = physSizeX*stackSizeX/newSizeX;
        newphysY = physSizeY*stackSizeY/newSizeY;
        newphysZ = physSizeZ;
        plane = uint16(zeros(newSizeX,newSizeY,newSizeZ,numchannels));
        
        
        if numchannels ==5
            for j = 1:stackSizeZ
                plane(:,:,j,1) = imresize(data{1,1}{(j*5)-4,1}, [newSizeX newSizeY]);
                plane(:,:,j,2) = imresize(data{1,1}{(j*5)-3,1}, [newSizeX newSizeY]);
                plane(:,:,j,3) = imresize(data{1,1}{(j*5)-2,1}, [newSizeX newSizeY]);
                plane(:,:,j,4) = imresize(data{1,1}{(j*5)-1,1}, [newSizeX newSizeY]);
                plane(:,:,j,5) = imresize(data{1,1}{(j*5)-0,1}, [newSizeX newSizeY]);
            end
        elseif numchannels == 4
            for j = 1:stackSizeZ
                plane(:,:,j,1) = imresize(data{1,1}{(j*4)-3,1}, [newSizeX newSizeY]);
                plane(:,:,j,2) = imresize(data{1,1}{(j*4)-2,1}, [newSizeX newSizeY]);
                plane(:,:,j,3) = imresize(data{1,1}{(j*4)-1,1}, [newSizeX newSizeY]);
                plane(:,:,j,4) = imresize(data{1,1}{(j*4)-0,1}, [newSizeX newSizeY]);
            end
        elseif numchannels == 3
            for j = 1:stackSizeZ
                plane(:,:,j,1) = imresize(data{1,1}{(j*3)-2,1}, [newSizeX newSizeY]);
                plane(:,:,j,2) = imresize(data{1,1}{(j*3)-1,1}, [newSizeX newSizeY]);
                plane(:,:,j,3) = imresize(data{1,1}{(j*3)-0,1}, [newSizeX newSizeY]);
            end
            
        elseif numchannels == 2
            for j = 1:stackSizeZ
                plane(:,:,j,1) = imresize(data{1,1}{(j*2)-1,1}, [newSizeX newSizeY]);
                plane(:,:,j,2) = imresize(data{1,1}{(j*2)-0,1}, [newSizeX newSizeY]);
            end
            
            
        elseif numchannels == 1
            for j = 1:stackSizeZ
                plane(:,:,j,1) = imresize(data{1,1}{(j*1)-0,1}, [newSizeX newSizeY]);
            end
            
        elseif numchannels == 8
            for j = 1:stackSizeZ
                plane(:,:,j,1) = max(imresize(data{1,1}{(j*8)-7,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*8)-3,1}, [newSizeX newSizeY]));
                plane(:,:,j,2) = max(imresize(data{1,1}{(j*8)-6,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*8)-2,1}, [newSizeX newSizeY]));
                plane(:,:,j,3) = max(imresize(data{1,1}{(j*8)-5,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*8)-1,1}, [newSizeX newSizeY]));
                plane(:,:,j,4) = max(imresize(data{1,1}{(j*8)-4,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*8)-0,1}, [newSizeX newSizeY]));
                
            end
            plane = plane(:,:,:,1:4);
            numchannels = 4;
            
        elseif numchannels == 10
            for j = 1:stackSizeZ
                plane(:,:,j,1) = max(imresize(data{1,1}{(j*10)-9,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*10)-4,1}, [newSizeX newSizeY]));
                plane(:,:,j,2) = max(imresize(data{1,1}{(j*10)-8,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*10)-3,1}, [newSizeX newSizeY]));
                plane(:,:,j,3) = max(imresize(data{1,1}{(j*10)-7,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*10)-2,1}, [newSizeX newSizeY]));
                plane(:,:,j,4) = max(imresize(data{1,1}{(j*10)-6,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*10)-1,1}, [newSizeX newSizeY]));
                plane(:,:,j,5) = max(imresize(data{1,1}{(j*10)-5,1}, [newSizeX newSizeY]),imresize(data{1,1}{(j*10)-0,1}, [newSizeX newSizeY]));
            end
            plane = plane(:,:,:,1:5);
            numchannels = 5;
        end
        %         metadata = createMinimalOMEXMLMetadata(uint16(zeros(size(plane))));
        
        names = {'pixel size X'; 'pixel size Y'; 'pixel size Z'; ...
            '# of Pixels X'; '# of Pixels Y'; '# of Pixels Z'; ...
            'old pixel size'; 'Total Channels'; 'Vessel Channel'; ...
            'Particle Channel'; 'Zoom'; 'LightSheet Thickness'; ...
            'Objective'; 'Laser 1 nm'; 'Laser 2 nm'; 'Laser 3 nm';...
            'Laser 4 nm'; 'Laser 5 nm'; 'Laser Power 1';...
            'Laser Power 2'; 'Laser Power 3';...
            'Laser Power 4'; 'Laser Power 5'; 'Filters 1'; 'Filters 2';...
            'Filters 3'; 'Filters 4'; 'Filters 5'; 'Exposure 1'; ...
            'Exposure 2'; 'Exposure 3'; 'Exposure 4'; 'Exposure 5';...
            'Track1'; 'Track2'; 'Track3'; 'Track4'; 'Track5'; 'ZLeft.405';...
            'ZLeft.488'; 'ZLeft.561'; 'ZLeft.638'; 'ZRight.405'; ...
            'ZRight.488'; 'ZRight.561'; 'ZRight.638' };
        newphys = {newphysX; newphysY; newphysZ; newSizeX; newSizeY; ...
            newSizeZ; physSizeX; numchannels; 3; numchannels; ...
            imagezoom; LSthickness; objective{1}; laser_wavelength(1);...
            laser_wavelength(2); laser_wavelength(3); laser_wavelength(4);...
            laser_wavelength(5); laser_powers(1); laser_powers(2); ...
            laser_powers(3); laser_powers(4); laser_powers(5);...
            filterset{1}; filterset{2}; filterset{3};...
            filterset{4}; filterset{5}; cam_int_time(1); cam_int_time(2);...
            cam_int_time(3);cam_int_time(4);cam_int_time(5); ...
            trackconfig{1};trackconfig{2};trackconfig{3};trackconfig{4};...
            trackconfig{5};leftrightalignment(1);leftrightalignment(2);...
            leftrightalignment(3);leftrightalignment(4);...
            leftrightalignment(5);leftrightalignment(6);...
            leftrightalignment(7);leftrightalignment(8);};
        
        T = table(names, newphys );
        writetable(T,metadatafilename);
        writetable(T2,metadatafullfilename);
        tic
        
        
        for k = 1:numchannels
            clear options
        options.overwrite = true;
        options.compress = 'no';
        saveastiff(plane(:,:,:,k),filenames(k,:),options)
        
%         imwrite(plane(:,:,1,k),filenames(k,:));
%             for p = 2:newSizeZ
%                 imwrite(plane(:,:,p,k),filenames(k,:), 'WriteMode','append');
%             end
%             
        end
        disp('writing time')
        toc
        
        %         save(metadatafilename,'stackSizeX','stackSizeY','stackSizeZ','numchannels','physSizeX','physSizeY','physSizeZ','newSizeX','newSizeY','newSizeZ','newphysX','newphysY','newphysZ')
    end
    cd(data_folder);
   % delete(files(int_file).name);
end

cd(matlab_folder)


% %%
function [laser_wavelength, laser_powers, cam_int_time, filterset, ... 
    trackconfig, leftrightalignment, LSthickness, imagezoom, objective]...
    = extractMetaData(metadatakeys_cell,metadatavalues_cell)
%%

        % Laser Wavelengths extraction
        index = [];
        for i = 1:5
            index = [index; find(strcmp(metadatakeys_cell, strcat('Global Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Attenuator|Wavelength #',num2str(i))))];
            
        end
        laser_wavelength= zeros(5,1);
        numlasers = size(index,1);
        for i= 1:numlasers
            laser_wavelength(i) = 1e9*str2num(metadatavalues_cell{index(i)});
        end

        % Laser Powers extraction
        index = [];
        for i = 1:5
            index = [index; find(strcmp(metadatakeys_cell, strcat('Global Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Attenuator|Transmission #',num2str(i))))];
            
        end
        laser_powers= zeros(5,1);
        for i= 1:size(index,1)
            laser_powers(i) = 100*str2num(metadatavalues_cell{index(i)});
        end
        
        % Camera Integration time extraction
        index = [];
        for i = 1:5
            index = [index; find(strcmp(metadatakeys_cell, strcat('Global Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|CameraIntegrationTime #',num2str(i))))];
            
        end
        cam_int_time= zeros(5,1);
        for i= 1:size(index,1)
            cam_int_time(i) = 1000*str2num(metadatavalues_cell{index(i)});
        end
        
        % Filter set extraction
        index = [];
        for i = 1:5
            index = [index; find(strcmp(metadatakeys_cell, strcat('Global Information|Instrument|Filter|Name #',num2str(i))))];
            
        end
        filterset= cell(5,1);
        for i= 1:size(index,1)
            filterset(i) = metadatavalues_cell(index(i));
        end
        
        %Track Configuration
        index = [];
        for i = 1:5
            index = [index; find(strcmp(metadatakeys_cell, strcat('Global Information|Instrument|Filter|Id #',num2str(i))))];
            
        end
        trackconfig= cell(5,1);
        for i= 1:size(index,1)
            trackconfig(i) = metadatavalues_cell(index(i));
        end

        %Left Right Configuration
        index = zeros(8,1);
        leftrightsuffixes = {'ZLeft.405'; 'ZLeft.488'; 'ZLeft.561'; 'ZLeft.638'; 'ZRight.405'; 'ZRight.488'; 'ZRight.561'; 'ZRight.638'};
        for i = 1:8
            index(i) = find(strcmp(metadatavalues_cell, strcat('CarlZeissAim.LightSheetOffset',leftrightsuffixes{i})));
            
        end
        index2 = [];
         for i = 1:8
            tempstr = metadatakeys_cell{index(i)};
            index2 = [index2;find(cell2mat(cellfun(@(x) size(regexp(char(x), strcat('\w*LsmTag #\w*',tempstr(end))),1), metadatakeys_cell, 'UniformOutput',false)),1)];
            
        end
        
        leftrightalignment= zeros(8,1);
        for i= 1:size(index2,1)
            leftrightalignment(i) = 100*str2num(metadatavalues_cell{index2(i)});
        end
        
        %Lightsheet thickness
        index1 = strcmp(metadatavalues_cell,'LightSheetThickness');
        tempstr = metadatakeys_cell{index1};
%         index2 = strcmp(metadatakeys_cell, strcat('Global LsmTag #',tempstr(end)));
        index2=find(cell2mat(cellfun(@(x) size(regexp(char(x), strcat('\w*LsmTag #\w*',tempstr(end))),1), metadatakeys_cell, 'UniformOutput',false)),1);
        LSthickness = str2num(metadatavalues_cell{index2});
        
        %Image zoom
        imagezoom = str2num(metadatavalues_cell{strcmp(metadatakeys_cell, 'Global Experiment|AcquisitionBlock|AcquisitionModeSetup|RtZoom #1')});        
        
        %Objective
        objective = metadatavalues_cell(strcmp(metadatakeys_cell, 'Global Information|Instrument|Objective|Manufacturer|Model #1'));

end