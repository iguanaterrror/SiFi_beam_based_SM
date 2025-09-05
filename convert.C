#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>

namespace fs = std::filesystem;

/**
 * @brief Converts a single ROOT file to a TSV file.
 * 
 * This function reads the "Secondaries" TTree from a ROOT file, filters
 * gamma particles (ParticleID == 22), applies a spatial cut, and writes
 * the particle information (ID, energy, position, momentum) to a TSV file.
 * 
 * @param input_path  Path to the input ROOT file.
 * @param output_path Path to the output TSV file.
 */
void root_to_tsv(const char* input_path, const char* output_path) {
    TFile* file = TFile::Open(input_path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "❌ Error opening file: " << input_path << std::endl;
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("Secondaries"));
    if (!tree) {
        std::cerr << "❌ Error: TTree 'Secondaries' not found in file.\n";
        file->Close();
        return;
    }

    // Pointers to store branch data
    std::vector<int>* ParticleID = nullptr;
    std::vector<double>* ParticleEnergy = nullptr;
    std::vector<TVector3>* ParticlePosition = nullptr;
    std::vector<TVector3>* ParticleMomentum = nullptr;
    std::vector<double>* ParticleTime = nullptr;

    // Link tree branches to pointers
    tree->SetBranchAddress("ParticleID", &ParticleID);
    tree->SetBranchAddress("ParticleEnergy", &ParticleEnergy);
    tree->SetBranchAddress("OriginPosition", &ParticlePosition);
    tree->SetBranchAddress("OriginMomentum", &ParticleMomentum);
    tree->SetBranchAddress("ParticleTime", &ParticleTime);

    Long64_t nentries = tree->GetEntries();
    std::ofstream outfile(output_path);
    if (!outfile.is_open()) {
        std::cerr << "❌ Error opening output file: " << output_path << std::endl;
        file->Close();
        return;
    }

    // Loop over all entries in the TTree
    for (Long64_t i = 0; i < nentries; ++i) {
        if (tree->GetEntry(i) <= 0) continue;

        // Only process gamma particles
        if ((*ParticleID)[0] == 22) {
            double x = (*ParticlePosition)[0].X();
            double y = (*ParticlePosition)[0].Y();
            double z = (*ParticlePosition)[0].Z();
            
            // Apply spatial cut (detector bounds)
            if (x >= -25 && x <= 25 && y >= -25 && y <= 25 && z >= -45 && z <= 45) {
                outfile << (*ParticleID)[0] << "\t"
                        << (*ParticleEnergy)[0] << "\t"
                        << x << "\t" << y << "\t" << z << "\t"
                        << (*ParticleMomentum)[0].X() << "\t"
                        << (*ParticleMomentum)[0].Y() << "\t"
                        << (*ParticleMomentum)[0].Z() << "\n";
            }
        }
    }

    outfile.close();
    file->Close();
}

/**
 * @brief Converts all ROOT files in subdirectories corresponding to energies.
 * 
 * Expects folders named like "Box_117dot8MeV_geom_containment". Extracts
 * the energy label from the folder name and writes TSV files to a corresponding
 * output directory.
 * 
 * @param main_input_dir  Path to the main input directory containing energy folders.
 * @param main_output_dir Path to the main output directory for TSV files.
 */
void convert_all_energies(const std::string& main_input_dir, const std::string& main_output_dir) {
    for (const auto& subdir : fs::directory_iterator(main_input_dir)) {
        if (!subdir.is_directory()) continue;

        std::string energy_folder = subdir.path().filename().string();
        if (energy_folder.rfind("Box_", 0) != 0) {
            std::cout << "⏭️ Skipping folder (does not start with 'Box_'): " << energy_folder << "\n";
            continue;
        }

        // Extract energy label from folder name
        size_t underscore_pos = energy_folder.find('_');
        size_t mev_pos = energy_folder.find("MeV");
        if (underscore_pos == std::string::npos || mev_pos == std::string::npos) {
            std::cerr << "⚠️ Unexpected folder name format: " << energy_folder << "\n";
            continue;
        }
        std::string energy_label = energy_folder.substr(underscore_pos + 1, mev_pos - underscore_pos - 1);

        // Create output directory for this energy
        std::string out_dir = main_output_dir + "/" + energy_label + "_extracted";
        fs::create_directories(out_dir);

        // Convert all ROOT files in this subdirectory
        for (const auto& file : fs::directory_iterator(subdir)) {
            if (file.path().extension() == ".root") {
                fs::path input_path = file.path();
                std::string original_stem = input_path.stem().string();
                fs::path output_path = fs::path(out_dir) / (original_stem + ".tsv");

                root_to_tsv(input_path.string().c_str(), output_path.string().c_str());
            }
        }
    }

    std::cout << "✅ All energies processed.\n";
}

/**
 * @brief Converts all ROOT files in a single folder.
 * 
 * This is useful when the folder does not follow the "Box_<energy>" naming scheme.
 * The folder name is used as the output label.
 * 
 * @param energy_folder_path Path to the folder containing ROOT files.
 * @param main_output_dir    Path to the main output directory for TSV files.
 */
void convert_one_folder(const std::string& energy_folder_path, const std::string& main_output_dir) {
    fs::path subdir(energy_folder_path);

    if (!fs::is_directory(subdir)) {
        std::cerr << "❌ Provided path is not a directory: " << energy_folder_path << "\n";
        return;
    }

    // Use folder name as energy label
    std::string energy_folder = subdir.filename().string();
    std::string energy_label = energy_folder;

    std::string out_dir = main_output_dir + "/" + energy_label + "_extracted";
    fs::create_directories(out_dir);

    for (const auto& file : fs::directory_iterator(subdir)) {
        if (file.path().extension() == ".root") {
            fs::path input_path = file.path();
            std::string original_stem = input_path.stem().string();
            fs::path output_path = fs::path(out_dir) / (original_stem + ".tsv");

            root_to_tsv(input_path.string().c_str(), output_path.string().c_str());
        }
    }

    std::cout << "✅ Folder processed: " << energy_folder << "\n";
}

/**
 * @brief Converts all subfolders in an input directory that match an optional prefix.
 * 
 * This is a flexible wrapper around convert_one_folder() for batch processing.
 * 
 * @param inputDir  Path to main input directory.
 * @param outputDir Path to main output directory.
 * @param prefix    Optional prefix to filter folders (default = "").
 */
void convert(const std::string& inputDir, const std::string& outputDir, const std::string& prefix = "") {
    for (const auto& subdir : fs::directory_iterator(inputDir)) {
        if (!subdir.is_directory()) continue;

        std::string folderName = subdir.path().filename().string();

        // Skip folders not matching the prefix
        if (!prefix.empty() && folderName.rfind(prefix, 0) != 0) {
            std::cout << "⏭️ Skipping folder (does not match prefix '" << prefix << "'): " << folderName << "\n";
            continue;
        }

        // Process this folder
        convert_one_folder(subdir.path().string(), outputDir);
    }
}
