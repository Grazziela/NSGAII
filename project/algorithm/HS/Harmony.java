package project.algorithm.HS;

import project.data.*;

import project.start.*;

import java.util.*;

public class Harmony 
{
	private double[] harmony = null;

	// The overall fitness of an individual (or harmony) is calculated based on how close an individual is to the desired sensitivity and specificity values, which are predefined by the user.
	private double sensitivity;
	private double specificy;
	private double fitness;

	// true positive, false negative, true positive and true negative variables used to calculate the sensitivity and specificity of the thresholds defined (true table)
	private int true_positive;
	private int false_negative;
	private int true_negative;
	private int false_positive;

	private int total;

	private int[] chromosome = null;

	// Class Harmony constructor method
	public Harmony(int length) 
	{
		// Harmony and chromosome have the same length. The chromosome has the information about the panel (binary). Harmony has the cut off values for screening.
		this.harmony = new double[length];
		this.chromosome = new int[length];
		initialise();
	}


	/*
	 * The objective of each harmony is to optimise the cut off for a certain Peptide/peptide during screening so that it optimises 
	 * the pair (Sensitivity,Specificity) for the whole data
	 * 
	 * In this version of the code, the harmony value is real and randomly selected from the cutoff windows, defined by the range 
	 * of values in the data set column.
	 */
	public void build(int[] chromosome, ArrayList<Peptide> peptides, RandomGenerator random_generator) 
	{
		this.chromosome = new int[chromosome.length];

		for (int i = 0; i < this.harmony.length; i++) // harmony and chromosome have the same length
		{
			this.chromosome[i] = chromosome[i];

			if (chromosome[i] == 1) // If the peptide is present in the panel
			{
				int index = random_generator.getRandomIndex(peptides.get(i).getCutOffsWindowSize());
				//if(index == 0) 
				//{
				//this.harmony[i] = 0.0;
				//} 
				//else 
				//{
				this.harmony[i] = peptides.get(i).getCutOffValue(index);
				//}
			}
		}
	}

	public void calculateVar(ArrayList<Sample> samples, int threshold) 
	{
		initialise();
		int cont;

		for (int i = 0; i < samples.size(); i++) 
		{
			Sample sample = new Sample();
			sample = samples.get(i);

			int count = 0;
			int panel_size = 0;

			for (int j = 0; j < this.harmony.length; j++) 
			{

				if (this.chromosome[j] == 1) // Peptide included in the panel
				{
					panel_size ++;

					if (sample.getPeptideValue(j) > this.harmony[j]) // Value above the cut off
					{
						if (sample.is_a_test_sample()) // if sample above the cut off is test, then increment true positives
						{
							this.true_positive++;
						} 
						else // increment false positives
						{
							this.false_positive++;
						}
						break;
					}

					else // value below cut off
					{
						count++;
					}
				}
			}

			if (count == panel_size) // all values for that sample are below the cut offs for the panel
			{
				if (!sample.is_a_test_sample()) // if the sample is control
				{
					this.true_negative++;
				} 
				else 
				{
					this.false_negative++;
				}
			}
		}
		if (this.true_positive + this.false_negative == 0) // this is to avoid division by zero
		{
			this.sensitivity = 0.0;
		} 
		else 
		{
			this.sensitivity = (double) this.true_positive/(double)(this.true_positive + this.false_negative);	
		}

		if (this.true_negative + this.false_positive == 0) 
		{
			this.specificy = 0.0;
		} 
		else 
		{
			this.specificy = (double)this.true_negative /(double)(this.true_negative + this.false_positive);
		}
	}

	public void calculateFitness(double sensitivity_goal, double specificity_goal, int generation, int total_generations) 
	{
		// Fitness calculation using the Chebchev method

		double scale = 1;
		double value = (generation+1)/scale;

		//double sens_weight = value/total_generations;
		//double spec_weight = (1 - sens_weight)/1.5;

		double sens_weight = Math.abs((sensitivity_goal - this.sensitivity)) * 1;
		double spec_weight = Math.abs((specificity_goal - this.specificy)) * 1.6;

		this.fitness = Math.max(Math.abs((sensitivity_goal - this.sensitivity)*sens_weight), Math.abs((specificity_goal - this.specificy)*spec_weight));
		//System.out.println("\n Fitness: " + this.fitness);
		this.total = this.false_negative+this.false_positive+this.true_negative+this.true_positive;
	}

	public void initialise() 
	{
		this.true_positive = 0;
		this.false_negative = 0;
		this.true_negative = 0;
		this.false_positive = 0;
		this.total = 0;
	}

	public void PrintOnScreen() 
	{
		if(this.true_positive >= 10 && this.true_negative >= 10) {
			System.out.println("The harmony value is " + Arrays.toString(this.harmony));
			System.out.println("The sensitivity is " + this.sensitivity + " and the specificity is " + this.specificy + " and the fitness value is " + this.fitness);
			System.out.println("The true positivie is: "+this.true_positive + "The true negative is: "+this.true_negative + "The false positive is: "+this.false_positive + "The false negative is: "+this.false_negative);
			System.out.println("The total number of Samples is " + this.total);
		}
	}

	public void setHarmony(Harmony harmony, int[] chromosome) 
	{
		setHarmonyValue(harmony.getHarmonyValue());
		setSensitivity(harmony.getSensitivity());
		setSpecificy(harmony.getSpecificy());
		setFitness(harmony.getFitness());
		setTruePositive(harmony.getTruePositive());
		setFalseNegative(harmony.getFalseNegative());
		setTrueNegative(harmony.getTrueNegative());
		setFalsePositive(harmony.getFalsePositive());
		setChromosome(chromosome);
	}

	public void setHarmonyIndexValue(int index, double value) 
	{
		this.harmony[index] = value;
	}

	public double getHarmonyIndexValue(int index) 
	{
		return this.harmony[index];
	}

	public void setChromosomeIndexValue(int index, int value) 
	{
		this.chromosome[index] = value;
	}

	public void setChromosome(int[] chromosome) 
	{
		for (int i = 0; i < chromosome.length; i++) 
		{
			this.chromosome[i] = chromosome[i];
		}
	}

	public double[] getHarmonyValue() 
	{
		double[] cloned = new double[this.harmony.length];
		for (int i = 0; i < cloned.length; i++) 
		{
			cloned[i] = this.harmony[i];
		}
		return cloned;
	}

	public void setHarmonyValue(double[] value) 
	{
		for (int i = 0; i < value.length; i++) 
		{
			this.harmony[i] = value[i];
		}
	}

	public double getSensitivity() 
	{
		return this.sensitivity;
	}

	public void setSensitivity(double sensitivity) 
	{
		this.sensitivity = sensitivity;
	}

	public double getSpecificy() 
	{
		return this.specificy;
	}

	public void setSpecificy(double specificy) 
	{
		this.specificy = specificy;
	}

	public double getFitness() 
	{
		return this.fitness;
	}

	public void setFitness(double fitness) 
	{
		this.fitness = fitness;
	}

	public void setTruePositive(int truePositive) 
	{
		this.true_positive = truePositive;
	}

	public int getTruePositive() 
	{
		return this.true_positive;
	}

	public void setTrueNegative(int trueNegative) 
	{
		this.true_negative = trueNegative;
	}

	public int getTrueNegative() 
	{
		return this.true_negative;
	}

	public void setFalsePositive(int falsePositive) 
	{
		this.false_positive = falsePositive;
	}

	public int getFalsePositive() 
	{
		return this.false_positive;
	}

	public void setFalseNegative(int false_negative) 
	{
		this.false_negative = false_negative;
	}

	public int getFalseNegative() 
	{
		return this.false_negative;
	}

	public void clear() 
	{
		this.harmony = new double[this.harmony.length];
	}

	public static Comparator<Harmony> sortByFITNESS = new Comparator<Harmony>() 
	{

		public int compare(Harmony harmony1, Harmony harmony2) 
		{

			double fitness_harmony1 = harmony1.getFitness();
			double fitness_harmony2 = harmony2.getFitness();

			/*For ascending order*/ 
			if (fitness_harmony1 - fitness_harmony2 > 0)
			{
				return 1;
			} 

			if (fitness_harmony1 - fitness_harmony2 < 0)
			{
				return -1;
			}

			return 0; 
		}
	};

	@Override
	public Harmony clone()
	{
		Harmony cloned = null;
		try 
		{
			cloned = (Harmony) super.clone();
		} catch (CloneNotSupportedException e) 
		{

		}
		return cloned;
	}
}
