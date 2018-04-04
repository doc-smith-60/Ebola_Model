/********************************************************************************
* Ebola_A.c																		*
* Created by Neal Smith on March 22, 2018.										*
* An attempt to streamline the processes from codes Ebola1.c through Ebola7.c.	*
* Code in this section should appear in working form, in the order of use.		*
*							Q:: Question.										*
*							F:: Needs fixing.									*
********************************************************************************/

//preprocessor directives
#include <stdio.h>	//Standard C input/output functions
#include <string.h>	//Library for functions on strings
#include <stdlib.h>	//Standard C Library (for memory allocation in particular) 
#include "Date_And_Reading_Reports.h"	//Header file for conversion of date to date_ID
#include "MTrandom.h"					//For random number generation (accept/reject situations)
#include "lfunc.h"						//For MTrandom.cpp

//definitions
#define MAX_REPORTS 120		//There are ~58000 entries in our dataset - only 111 in the Mali subset.
								//Likely value: 100000. Reduce this number if memory issues occur
#define NUMDAYS 428


//global constants, variables and structures
	//Anything beginning with "p_" is a pointer.
char code_name[100];
char case_file_name[100];

struct parameter_list *p_parameters = &parameters;		//A pointer initialised to a parameter list
struct current_case_report *p_reports = &report_list;	//A pointer initiliased to a case report list
	//POTENTIAL ISSUE: REPORT LIST IS AN ARRAY OF STRUCTURES. I'M NOT SURE HOW BEST TO REFERENCE IT.

	//Structures defining the following located in "Date_And_Reading_Reports.h":
		//An individual case
		//A case report
		//Model parameters
	//This file really needs a new name

//functions
/*-----------------------------------
|FUNCTIONS REGARDING ACCESS TO DATA |
-----------------------------------*/
//2) Warning if anything is wrong with files passed to function
void usage()
{
	printf("Command line should contain the following files, each with names < 100ch:\nCase data input file.\n");
}

//1) To ensure we have the files we need.
void handleargs(int argc, char **argv)
{
	printf("Determining if all necessary files are present.\n");
	if (argc != 2) usage();													//Number of files in command line, plus one(for the filename). We currently have one file - the case data.
	else printf("Correct number of files provided as command arguments.\n");
	if (strlen(argv[1])>100) usage();										//Making sure the title isn't too long - remnant of Jon's code. //Q:: needed?
	strcpy(case_file_name, argv[1]);										//Copy the file name for a check, and to call file from.
	printf("Case file name is %s.\n", case_file_name);							//To ensure correct files are in correct locations.
}

/*--------------------------------------------------
| FUNCTIONS REGARDING READING OF DATA INTO PROGRAM |
--------------------------------------------------*/

//3) Read case data into an array, so it can be used to generate a linked list of cases
int read_case_data(struct parameter_list *p_params)	//F::check naming convention for called/calling fns
{
	//local variables
	FILE *patient_data;		//For the input file
	char line[100];			//To read unnecessary entries to (e.g. column headings)
	char ch;				//To read individual characters of names to temporarily
	int i;					//To read country names into the report list
	int j;					//To read subregion names into the report list
	int n;					//For the line in the report list that we are up to
//	int m;					//To count the number of cases that have been created
	
	printf("Opened read_case_data.\n");

	//Allocate sufficient memory to store the number of case reports
	report_list = (struct current_case_report*)malloc(MAX_REPORTS * sizeof(struct current_case_report));
	if (!report_list){
		printf("Could not allocate report_list.\n");
		exit(1);
	}
	else {
		printf("Allocated report_list successfully.\n");
	}

	//Use data in file on command line and report_list to store patients
	printf("Reading case data into array.\n");

	patient_data = fopen(case_file_name, "r");
	if (patient_data == NULL)
	{
		printf("\n File containing patient data could not be opened\n");
		exit(1);
	}
	else printf("Patient data file is open.\n");
	fgets(line, 1000, patient_data);
				//reads header. Ready to read first line of case data.
	
	//Setting up variables to read data
	p_params->num_diagnosed = 0;			//Reported case count
	p_params->total_cases = 0;				//Total cases
	n = 0;									//First line of report list
	
	//Reading data into array
	do
	{
		//Read country into array
		i = 0;
		report_list[n].country[i] = 0;
		while (fscanf(patient_data, "%c", &ch) == 1 && ch != ',' && i < 50)	//i<50 to match maximum length of country name
			report_list[n].country[i++] = ch;								//Read up to and including first comma (country name)
		report_list[n].country[i] = 0;										//Finish string for country name
		//Read subregion into array
		j = 0;
		while (fscanf(patient_data, "%c", &ch) == 1 && ch != ',' && j < 100)//j<100 for same reason as i<50 above
			report_list[n].subregion[j++] = ch;								//Read up to and including second comma (subregion name)
		report_list[n].subregion[j] = 0;									//Finish string for subregion name

		//Read number of cases in report
		fscanf(patient_data, "%d", &report_list[n].cases);
		fscanf(patient_data, "%c", &ch);		//reads third comma

		//Read date of report
		fscanf(patient_data, "%d/%d/%d", &report_list[n].day, &report_list[n].month, &report_list[n].year);
		report_list[n].date_ID = generate_report_date_id(report_list[n].day, report_list[n].month, report_list[n].year);

		//If there are cases on this date, generate the cases in a linked list
			if (report_list[n].cases >= 1){
			for(i = 1; i <= report_list[n].cases; i++){		//For all cases in this report,
				generate_cases(p_reports, &head, p_params);	//Make a case
			}	//Check to see if passing these are needed
					//POTENTIAL ISSUE: MIXING UP OF USING p_reports AND report_list
						//report_list is a pointer to the start of an array. Can I pass it to a fn?
							//If not, should I pass p_reports to this (and all) fns instead?
		printf("Generating case for next report.\nGenerated cases = %d.\n", p_params->num_diagnosed);	
		}

		//Check reading of entry and determine end of reports
		if (report_list[n].country[0] != 0) {	//If there is data in the country column (i.e. there is an entry)
			printf("Data for report %d:\n\tCountry: %s\n", n + 1, report_list[n].country);
			printf("\tSubregion: %s\n", report_list[n].subregion);
			printf("\tCases: %d\n", report_list[n].cases);
			printf("\tReport_Date: %d/%d/%d\n", report_list[n].day, report_list[n].month, report_list[n].year);
			printf("\tReport_Date_ID: %d\n", report_list[n].date_ID);
			n++;									//Ready to enter next report in next array 'row'
			p_params->total_reports = n;		//Update number of reports
		}
		else {	//If there is no entry,
			p_params->total_reports = n;
			printf("Reading of reports complete.\nTotal number of reports = %d.", p_params->total_reports);
		}
	} while (fgets(line, 1000, patient_data));			//Reads rest of line, while there are lines
	fclose(patient_data);			//close file when data read from it

	return 0;
}
/*---------------
| MAIN FUNCTION |
---------------*/

int main(int argc, char **argv)
{
	strcpy(code_name, argv[0]);
	printf("Starting code for %s.", code_name);		//Baseline code - so something is happening!
	
	//Make sure the necessary input files are present
	handleargs(argc, argv);

	//Read case report data into an array
		//Step two: convert data into cases. Not yet.
	read_case_data(p_parameters);

	//Initialise parameters

	//Open output file and set up counters for MCMC

	//Initiate Gibbs sampler

	getchar();					//So I can see what I've done - otherwise the code exits

	return 0;					//Gives an integer to the computer so the function has returned a value, as indicated.
}