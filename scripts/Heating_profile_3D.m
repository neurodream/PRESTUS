%% Set-up
clear;
heating_profile = niftiread("D:\heating_results\sub-003_final_heating_active2.nii.gz");

%% find max value
max_val = max(heating_profile(:));
linear_ind = find(heating_profile == max_val, 1);
[x,y,z] = ind2sub(size(heating_profile), linear_ind);

% max value is in skull.

%%
volshow(heating_profile)

%% Adjusting heating profile
heating_profile_adjusted = heating_profile - 37; % New variable to prevent further decreasing when re-running section

min_temp = min(heating_profile_adjusted, [], "all"); % check min
sprintf("Minimum after substracting 37: %d", min_temp)

figure; % Create a new figure
histogram(heating_profile_adjusted); %Visualize data to check for outliers
title('Visualizing Deviations from 0');
xlabel('Value');
ylabel('Frequency');

%%
volshow(heating_profile_adjusted)

%% Removing Outliers
heating_profile_adjusted(heating_profile_adjusted < 6) = 0;

% removing the most prominent outliers very conservatively
min_temp = min(heating_profile_adjusted, [], "all");
sprintf('Min after removing outliers: %d', min_temp)

figure; % Create a new figure
histogram(heating_profile_adjusted); %Plot again
title('Visualizing Deviations from 0');
xlabel('Value');
ylabel('Frequency');

%%
volshow(heating_profile_adjusted, 'Colormap', hot(256))

%% Assumption of 0 at voxels without heating
heating_profile_adjusted(heating_profile_adjusted < 0) = 0;
%Apparently there are more outliers with out a visible peak in the historam ... Thus, less conservative corrections

min_temp = min(heating_profile_adjusted, [], "all");
sprintf('Sanity check: This should be 0 -->: %d', min_temp)

figure; % Create a new figure
histogram(heating_profile_adjusted); %Plot again
title('Visualizing Deviations from 0');
xlabel('Value');
ylabel('Frequency');

%%
volshow(heating_profile_adjusted, 'Colormap', hot(256))

% The assumption seems to make sense! 

%%
volshow(heating_profile, 'Colormap', hot(256), 'VolumeThreshold', 0)


%% Introducing the MRI Image

% Subject003: T1-weighted
mri_scan = niftiread("/Users/lutztebbe/Documents/Master_Internship/TUS_Simulations/Subject_data/scans/sub-003_T1.nii.gz");

volshow(mri_scan)


