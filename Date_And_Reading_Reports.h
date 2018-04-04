/********************************************************************************
*	Date_And_Reading_Reports.h													*
*	Created by Neal Smith in March 2018											*
*	Contains:																	*
*		- Structures required for both Ebola_x.c and Date_And_Reading_Reports.c	*
*		- Functions defined in Date_And_Reading_Reports.c						*
********************************************************************************/

//Definitions of global variables
#define NUMDAYS 428			//The number of days from the index case being infectious to the last report
								//I need to verify this ASAP

/********************************************
* Structures required for both source files *
********************************************/

//Structure describing case
typedef struct patient* p_patient;		//Q:: Allowing p_patient to be the data type "pointer to patient"?
struct patient
{
	int index;		//Used to identify case position in array for updating founder
	char country[50];		//Used to generate the specific location
	char subregion[100];	//Used to generate the specific location
	double x;				//x co-ordinate - generated from the above
	double y;				//y co-ordinate - generated from the above
	int dates[4];		//key dates for a case:
						//dates[0]: day of exposure(or ?seeding event?)
						//dates[1]: day first infective(showing symptoms)
						//F:: dates[0] or dates[1] = -1: ?initial? case (possibly index)
						//dates[2]: date of death or symptom resolution
						//dates[3]: date of burial(or end of infectivity)
						//F:: We either need to segregate all three transmission types here and then write code linking these dates to the transmisison,
						//or we need to have an extra variable for the state a patient is in. I	personally prefer the first option.
						//specifics of these dates outlined in thesis
	p_patient parent_case;	//case of origin: NULL for cases with dates[0] = -1 or dates[1] = 1
								//This distinction depends on what we choose to be the start date
	int transmission_type;	//four types eventually, two to start with:
							//0 for index transmission (or zoonotic - just index case)
							//1 for individual transmission (to begin with, all non-index cases)
							//2 for hospital / care transmission
							//3 for burial transmission
	int survive;			//1 for survival, 0 otherwise
	int est_case;			//1 for seeding event resulting in case, 0 otherwise
	int diag;				//0 for no diagnosis (?unknown?	case), 1 for diagnosis
	int diag_day;			//date of diagnosis (9999 for unknown)
							//This will be the date chosen when creating a case
	int secondary_cases;	//from this case - can be computed but easier to store (from equivalent in Jon's code)
	int secondary_cases_gen[NUMDAYS]; //number of secondary cases exposed ('seeding events') each day by this case
	p_patient first_2dary;	//one case that this case exposed - so we can make a list of "children"
	p_patient left_sib;		//other case exposed by same primary case
	p_patient right_sib;	//other case exposed by same primary case
							//int search_type[NUMDAYS+1];	//Stores the search type for each generation - MAY BE USEFUL FOR INTERVENTIONS ?
	int pop_dens;			//Stores the population density for the case location (we may assume constant)
							//int treated[NUMDAYS+1];		//Stores whether subject to treatment - MAY BE USEFUL FOR INTERVENTIONS ?
							//int num_unobs ;				//Number of unobserved 2dary cases - IF WE WANT TO STORE THIS
	p_patient next;
	p_patient prev;
} *head;					//Create the first case in the linked list

	//parameters of model
struct parameter_list
{
	int total_reports;		//Total number of case reports - USED TO CHECK READING OF CASE DATA
	int num_diagnosed;		//Should come from reading the case reports
	int total_cases;		//Included num_diagnosed and undiagnosed
} parameters;				//(pointer to?) parameter list used in model

	//Structure of the date from each case report
struct current_case_report		//Because we don't have a line list, read each report's data here first
{								//From there, create case list based on data smoothing techniques agreed upon.
	char country[50];		//Lengths of these strings longer than needed on purpose
	char subregion[100];		//Country and subregion for each report - BE CAREFUL OF NATIONAL REPORTS AS WELL
	char report_type[30];	//new cases, cumulative cases, deaths etc.
							//F:: NOT IN CURRENT DATASET. EXPAND UPON FOR LARGER DATASET.
	int cases;				//How many cases were on that day, in that region
	int day;				//Day in the date of report
	int month;				//Month in the date of report
	int year;				//Year in the date of report
	int date_ID;			//The date of this report.
	int previous_date_ID;	//The date of the previous report FOR THE SAME SUBREGION - used for temporal precision
} *report_list;	 //F:: Figure out how to best access these reports when needed.

/****************************************
* Functions defined in this source file *
****************************************/

int generate_report_date_id(int day, int month, int year);
void generate_cases(struct current_case_report *current_report, p_patient *head, struct parameter_list *p_params);
int assign_x();
int assign_y();
