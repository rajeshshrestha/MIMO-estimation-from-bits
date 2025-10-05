function []=main(method_type, data_type, N_trials, Nr, Nt, K_true,R, quant_bits,K_model, test_name, data_file_path)

% Main script to run the experiments
% method_type: 'admm', 'newadmm', 'all'
% data_type: 'random', 'from_file'
% N_trials: number of trials to run
% Nr: number of receiver antennas
% Nt: number of transmitter antennas
% K_true: number of paths in the true channel
% R: number of measurements to be sent back
% quant_bits: number of bits to quantize each measurement
% K_model: number of paths in the model used for recovery
% test_name: name of the test to save data
% data_file_path: path to the data file if data_type is 'from_file' [only required if data_type is 'from_file']

% Randomize
% rng('shuffle');
rng(42); % For reproducibility during testing


% Add necessary files to path
mfilepath = fileparts(which(mfilename));
addpath(fullfile(mfilepath, 'dataGenerators'));
addpath(fullfile(mfilepath, 'HFeedbacks'));
addpath(fullfile(mfilepath, 'methods'));
addpath(fullfile(mfilepath, 'utils'));


% Change the variables to char
method_type = convertCharsToStrings(method_type);
data_type = convertCharsToStrings(data_type);


% Define parameters
alpha = 0.5;

HSNR=str2double(getenv("HSNR"));
if isnan(HSNR)
    HSNR = 300; % SNR value in dB for quality of channel present in ME
end

BER = str2double(getenv("BER"));
if isnan(BER)
    BER = 0; % Bit error rate for the feedback channel
end 

YSNR = str2double(getenv("YSNR"));
if isnan(YSNR)
    YSNR = NaN; % SNR value in dB for quality of received pilot signals in ME
end 

Npilotsfactor = str2double(getenv("Npilotsfactor"));
if isnan(Npilotsfactor)
    Npilotsfactor = 3; % how many times the Nt should be for number of pilots
end 

SNR=str2double(getenv("DitherSNR"));
if isnan(SNR)
    SNR = 30; %SNR value in dB for dithering
end



if data_type == "random"
    data_generator = RandomHGenerator(K_true, Nr, Nt, alpha, HSNR);
elseif data_type == "from_file"
    data_generator = FileHGenerator(K_true, Nr, Nt, HSNR, data_file_path);
else
    error("Unknown data generator passed: "+data_type)
end

feedback_generator = FeedbackGenerator(SNR, quant_bits, BER, YSNR);


% Intialize lists to store actual channels
H_true_lst = [];
H_admm_lst = [];
H_new_admm_lst = [];


% Initialize lists to store estimated channels
NMSE_admm_lst =[];
NMSE_new_admm_lst = [];

% Intialize lists to store feedback bit numbers
feedback_bitnum_admm_lst =[];
feedback_bitnum_new_admm_lst =[];



% Turn off existing directory warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');
for n_trial=1:N_trials
    fprintf("-----------------------------------------------------------------------\n")
    fprintf("Running Trial: %d \n", n_trial)
    fprintf("-----------------------------------------------------------------------\n")

    % Generate channel matrix and its vectorised real h
    [H_true,H_noisy] = data_generator.getH();
    H_true_lst = cat(1, H_true_lst,reshape(H_true, [1,Nt,Nr]));

    if method_type == "admm" || method_type == "all"
        [z, A, sigma_true] = feedback_generator.feedbackADMM(H_true, R, Npilotsfactor);
        H_admm = admm(z,A,sigma_true,feedback_generator.scaler_for_quantization,Nr,Nt,K_model);
        NMSE_admm = norm(H_admm/norm(H_true,'fro') - H_true/norm(H_true,'fro'),'fro')^2;
        fprintf("NMSE for ADMM: %f \n", NMSE_admm)

        H_admm_lst = cat(1, H_admm_lst,reshape(H_admm, [1,Nt,Nr]));
        NMSE_admm_lst(end+1) = NMSE_admm;
        feedback_bitnum_admm_lst(end+1) = R * quant_bits;
    end
    if method_type == "newadmm" || method_type == "all" % New compression technique
        [z_new, A_new, sigma_true_new, ] = feedback_generator.feedbackNewADMM(H_true,R);
        H_new_admm = newADMM(z_new,A_new,sigma_true_new,feedback_generator.scaler_for_quantization,Nr,Nt,K_model);
        NMSE_new_admm = norm(H_new_admm/norm(H_true,'fro') - H_true/norm(H_true,'fro'),'fro')^2;
        fprintf("NMSE for New Compression ADMM: %f \n", NMSE_new_admm)

        H_new_admm_lst = cat(1, H_new_admm_lst,reshape(H_new_admm, [1,Nt,Nr]));
        NMSE_new_admm_lst(end+1) = NMSE_new_admm;
        feedback_bitnum_new_admm_lst(end+1) = R * quant_bits;
    end
    fprintf("=======================================================================\n\n")

end
% Data save parent folder
parent_folder = mfilepath;
data_save_dir = fullfile(parent_folder,"exp-data",test_name,"data-type="+data_type+"/Nt="+num2str(Nt)+"_Nr="+num2str(Nr)+"/Ktrue="+num2str(K_true)+"_K_model="+num2str(K_model));
mkdir(data_save_dir)

fprintf("***********************************************************************\n")
fprintf("Average NMSE\n")
fprintf("***********************************************************************\n")
% Generate file name
file_base_name = "R=" + num2str(R) + "_" + "quant_bits=" + num2str(quant_bits) + "_" + "SNR=" + num2str(HSNR)+"_" + num2str(randi(10e6,1))+ "_" + num2str(randi(10e6,1))+ "_" + num2str(randi(10e6,1));
if method_type == "admm" || method_type == "all"

    fprintf("Average NMSE for ADMM: %f \n", mean(NMSE_admm_lst))
    method_folder = fullfile(data_save_dir,"method=admm");
    mkdir(method_folder)

    file_num = numel(dir(fullfile(method_folder,"*.mat")))+1;
    file_name = file_base_name + "_" + num2str(file_num) + ".mat";
    file_path = fullfile(method_folder, file_name);
    save(file_path, "H_true_lst", "H_admm_lst","NMSE_admm_lst", "feedback_bitnum_admm_lst","method_type", "data_type", "N_trials", "Nr", "Nt", "K_true", "R", "quant_bits", "K_model", "alpha", "SNR","HSNR","BER", "YSNR")
end
if method_type == "newadmm" || method_type == "all"

    fprintf("Average NMSE for New Compression ADMM: %f \n", mean(NMSE_new_admm_lst))
    method_folder = fullfile(data_save_dir,"method=new_admm");
    mkdir(method_folder)

    file_num = numel(dir(fullfile(method_folder,"*.mat")))+1;
    file_name = file_base_name + "_" + num2str(file_num) + ".mat";
    file_path = fullfile(method_folder, file_name);
    save(file_path, "H_true_lst", "H_new_admm_lst","NMSE_new_admm_lst", "feedback_bitnum_new_admm_lst","method_type", "data_type", "N_trials", "Nr", "Nt", "K_true", "R", "quant_bits", "K_model", "alpha", "SNR","HSNR","BER", "YSNR")
end

fprintf("=======================================================================\n")
end