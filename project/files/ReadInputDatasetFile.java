package project.files;

import project.data.*;

import java.io.*;
import java.util.*;

public class ReadInputDatasetFile
{

	private ArrayList<Sample> Sample_Values = null; // this list stores the values for each peptide per sample
	private ArrayList<Peptide> Peptide_Names = null; // This list stores data file column headers containing the names of the peptides/antigens

	public ReadInputDatasetFile(String input_data_file_name) 
	{
		this.Sample_Values = new ArrayList<Sample>();
		this.Peptide_Names = new ArrayList<Peptide>();

		String information = null;
		int line_number = 0;

		// Reading input data file
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(input_data_file_name));

			while ((information = br.readLine()) != null) 
			{
				if (!information.trim().equals("")) 
				{
					String[] info = information.split(",");
					int number_Peptides = info.length;

					/* If the line is equal to zero, that means the program is reading the dataset headers. 
					 *  The headers are the names of the peptides/antigens.
					 *  These names need to be kept for when we decodify the chromossome and present the panels back to the researchers.
					 *  The number of peptides/antigens will be equal to the number of columns in the file minus one. 
					 *  The last column contains the class labels of the samples
					 */
					if (line_number == 0) 
					{
						// Capturing all names of the antigens/peptides in the dataset columns and saving in the list array.
						for (int i = 0; i < number_Peptides - 1; i++) 
						{
							Peptide peptide = new  Peptide(info[i+1], i);
							this.Peptide_Names.add(peptide);
						}
					} 
					/* Now reading the measured values for each peptide/antigen for a sample.
					 * The input file has the format: rows are samples, columns are the antigens/peptides, 
					 * last column is the class
					 */
					else 
					{
						// A sample has an id (sequential number), and the peptides' values
						Sample sample = new Sample(line_number - 1, info);					
						this.Sample_Values.add(sample);

						for (int i = 0; i < number_Peptides - 1; i++) 
						{
							this.Peptide_Names.get(i).update(Double.parseDouble(info[i]));
						}

					}
					line_number++;
				}
			}
			br.close();
		} 
		catch(IOException ioe) {
			System.err.println("Scenario File Reading Error!!!!!");
			System.exit(-1);
		}
	}

	public ArrayList<String[]> getFileInformation(String file_name) {
		String information = null;
		ArrayList<String[]> file_information = new ArrayList<String[]>();

		try{
			BufferedReader br = new BufferedReader(new FileReader(file_name));
			while ((information = br.readLine()) != null) {
				if (!information.trim().equals("")) {
					String[] info = information.split(",");
					file_information.add(info);
				}
			}

			br.close();

		} catch(IOException ioe) {
			System.err.println("Scenario File Reading Error!!!!!");
			System.exit(-1);
		}

		return file_information;
	}

	public ArrayList<Peptide> getPeptideNames()
	{
		// Generating the array with pre-defined set of feasible cut offs for Harmony Search
		for (int i = 0; i < this.Peptide_Names.size(); i++)
		{
			this.Peptide_Names.get(i).PreDefiningFeasibleCutOffs();
		}
		return this.Peptide_Names;
	}

	public ArrayList<Sample> getSampleValues()
	{
		return this.Sample_Values;
	}
	
	public int getChromosomeSize()
	{
		return this.Peptide_Names.size();
	}
}
