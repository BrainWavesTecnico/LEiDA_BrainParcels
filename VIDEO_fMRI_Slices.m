function VIDEO_fMRI_Slices

%%%%%%%%%%%%%%%%%%
%
%  Script to make a video of the fMRI signals over time
%
%
%%%%%%%%%%%%%%%%%%

% Choose scan, unzip and load to Workspace
% USER ADAPT TO THE EXTENSION WHERE YOU HAVE A .nii file
fMRI_file='/Users/joana/Documents/Work/Conferences/NeurimageCourse/SPMexercise/SPMexercise/swrrasub-MRI201906051BART007_task-reappraisal_run-01_bold.nii';
%gunzip(fMRI_file); %uncomment if it is a .nii.gz
fMRI=niftiread(fMRI_file);

% Check the size of the file:
[X, Y, Z, T]=size(fMRI); % Save size

% Define orientation of the slices
Orientation='Axial'; % or 'Sagital'

if strcmp(Orientation,'Axial')
    slices_to_plot=round((Z/11):(Z/11):Z-(Z/11)); %[5 10 15 20 25 30 35 40 45 50 55 60];
    fMRI=single(fMRI(:,:,slices_to_plot,:)); % Keep only the slices that we chose
    Z=length(slices_to_plot);
elseif strcmp(Orientation,'Sagital')
    slices_to_plot=round((X/11):(X/11):X-(X/11)); %[5 10 15 20 25 30 35 40 45 50 55 60];
    fMRI=single(fMRI(slices_to_plot,:,:,:)); % Keep only the slices that we chose
    X=length(slices_to_plot);
end

% Decide if you want to save video (debug first with save_video=0)
save_video=1;
name_video='fMRI_prepro_sagital'; % Define here the name the video will be saved

% Re-Organize the signals to have a size NxT (Voxels vs Time)
fMRI=detrend(reshape(fMRI,[X*Y*Z T]));
fMRI=fMRI-mean(fMRI,2); % Remove the mean
fMRI=reshape(fMRI,[X, Y, Z, T]); % Re-organize to 4D format

%%%% MAKE THE VIDEO

figure('color','k','Position',[75   377  1366   420])
colormap(jet)

if save_video
    videoModes = VideoWriter([name_video '.mp4'],'MPEG-4'); %create the video object
    videoModes.FrameRate = 5;  % number of frames per second
    videoModes.Quality = 100;
    open(videoModes); %open the file for writing
end

% Define a value for the thresold, to use the same in all images.
lim_colorbar=5*std(abs(fMRI(:)));

for t=1:T
    
    for slice = 1:length(slices_to_plot)
        
        subplot_tight(1,length(slices_to_plot),slice)
        
        if strcmp(Orientation,'Axial')
            imagesc(squeeze(fMRI(:,:,slice,t))',[-lim_colorbar lim_colorbar])
        elseif strcmp(Orientation,'Sagital')
            imagesc(squeeze(fMRI(slice,:,:,t))',[-lim_colorbar lim_colorbar])
        end
        
        if slice==1
            title({['TR=' num2str(t)],'',['Slice=' num2str(slices_to_plot(slice))]},'color','w','interpreter','none')
        else
            title(['Slice=' num2str(slices_to_plot(slice))],'color','w')
        end
        
        axis equal
        axis xy
        if strcmp(Orientation,'Axial')
            xlim([0 X])
            ylim([0 Y])
        elseif strcmp(Orientation,'Sagital')
            xlim([0 X])
            ylim([0 Z])
        end
        axis off
    end
    
    pause(0.1)    %  pause
    
    if save_video
        frame = getframe(gcf);
        writeVideo(videoModes,frame); %write the image to file
        clf
    end
    
end

if save_video
    close(videoModes); %close the file
end
