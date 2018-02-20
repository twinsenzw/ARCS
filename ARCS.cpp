#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <stdexcept>
#include <algorithm>
using namespace::std;

class genotype
{
public:
	vector<vector<string>> alleletypes; // chr x site x 2types
	vector<vector<int>> positions; // chr x site
	vector<vector<string>> sitenames;
};


struct location // location of a snp in the chromosome map, use as map.second
{
public:
	int chr;
	int bp;
	double cM;
};

double phi(double x) //cdf of standard normal
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

void ParseMap(ifstream &snpTraits, map<string, location> &snpLocations)
{
	string line;
	for(;getline(snpTraits,line);)
	{
		if(line=="") continue;
		int chr;
		string ID;
		double cM;
		int bp;
		
		istringstream linestream(line);
		linestream >> chr >> ID >> cM >> bp;
		location newlocation;
		newlocation.chr=chr;
		newlocation.cM=cM;
		newlocation.bp=bp;
		snpLocations[ID]=newlocation;
	}
}

void ParseInput(ifstream &input, map<string, location> snpLocations, map<string, genotype> &genotypes)
{
	string line;
	vector<string> sampleNames;
	vector<genotype> sampleGenotypes;

	getline(input,line);
	istringstream linestream1(line);
	string field;
	getline(linestream1,field,',');
	getline(linestream1,field,',');
	//getline(linestream1,field,',');

	for(;getline(linestream1,field,',');) sampleNames.push_back(field);

	sampleGenotypes.resize(sampleNames.size());
	for(auto sampleGenotype=sampleGenotypes.begin();sampleGenotype!=sampleGenotypes.end();sampleGenotype++)
	{
		sampleGenotype->alleletypes.resize(22);
		sampleGenotype->positions.resize(22);
		sampleGenotype->sitenames.resize(22);
	}
	for(;getline(input,line);)
	{
		if(line=="") continue;
		istringstream linestream(line);
		string snpID;
		getline(linestream,field,',');
		//getline(linestream,field,',');
		getline(linestream,snpID,',');

		auto snpSearch=snpLocations.find(snpID);
		
		if(snpSearch!=snpLocations.end())
		{
			//cout<<"check"<<endl;
			vector<string> allele;
			for(int samplei=0;getline(linestream,field,',');samplei++)
			{
				allele.push_back(field);
			}
			if(allele.size()!=sampleNames.size()) continue; //some samples don't have all records for specific sites

			for(int samplei=0;samplei<allele.size();samplei++)
			{
				
				sampleGenotypes[samplei].alleletypes[snpSearch->second.chr-1].push_back(allele[samplei]); //actually sorting the characters does not help much.
				sampleGenotypes[samplei].positions[snpSearch->second.chr-1].push_back(snpSearch->second.bp);
				sampleGenotypes[samplei].sitenames[snpSearch->second.chr-1].push_back(snpID);
			}
		}
				
	}

	for(int cnt=0;cnt<sampleNames.size();++cnt)
	{
		genotypes[sampleNames[cnt]]=sampleGenotypes[cnt];
	}
}

void comparison(string sample1, string sample2, int chr, map<string, genotype> genotypes, ofstream &output) //output genome comparison graph
{
	int chrindex=chr-1;
	for(int sitei=0;sitei<genotypes[sample1].positions[chrindex].size();sitei++)
	{
		int numMatchMax=0;
		string s1Alleletype=genotypes[sample1].alleletypes[chrindex][sitei];
		string s2Alleletype=genotypes[sample2].alleletypes[chrindex][sitei];

		if(s1Alleletype.find("-")!=string::npos || s2Alleletype.find("-")!=string::npos) continue;
		if(s1Alleletype.size()==0 || s1Alleletype.size()==0) continue;

		if(s1Alleletype[0]==s2Alleletype[0])
		{
			if(s1Alleletype[1]==s2Alleletype[1]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}	
	
		if(s1Alleletype[0]==s2Alleletype[1])
		{
			if(s1Alleletype[1]==s2Alleletype[0]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}

		output <<genotypes[sample1].positions[chrindex][sitei]<<"\t"<<numMatchMax<<endl;
	}
}

void comparison_extract0(string sample1, string sample2, int chr, map<string, genotype> genotypes, ofstream &output, vector<int> &zero) //extract sites that have zero matches
{
	int chrindex=chr-1;
	for(int sitei=0;sitei<genotypes[sample1].positions[chrindex].size();sitei++)
	{
		int numMatchMax=0;
		string s1Alleletype=genotypes[sample1].alleletypes[chrindex][sitei];
		string s2Alleletype=genotypes[sample2].alleletypes[chrindex][sitei];

		if(s1Alleletype.find("-")!=string::npos || s2Alleletype.find("-")!=string::npos) continue;

		if(s1Alleletype[0]==s2Alleletype[0])
		{
			if(s1Alleletype[1]==s2Alleletype[1]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}	
	
		if(s1Alleletype[0]==s2Alleletype[1])
		{
			if(s1Alleletype[1]==s2Alleletype[0]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}

		if(numMatchMax==0) zero.push_back(genotypes[sample1].positions[chrindex][sitei]);
	}
}

void comparison_usefilter(string sample1, string sample2, int chr, map<string, genotype> genotypes, ofstream &output, vector<int> filter)
{
	int chrindex=chr-1;
	for(int sitei=0;sitei<genotypes[sample1].positions[chrindex].size();sitei++)
	{
		if(find(filter.begin(),filter.end(),genotypes[sample1].positions[chrindex][sitei])==filter.end()) continue;
		int numMatchMax=0;
		string s1Alleletype=genotypes[sample1].alleletypes[chrindex][sitei];
		string s2Alleletype=genotypes[sample2].alleletypes[chrindex][sitei];


		if(s1Alleletype.find("-")!=string::npos || s2Alleletype.find("-")!=string::npos) continue;

		if(s1Alleletype[0]==s2Alleletype[0])
		{
			if(s1Alleletype[1]==s2Alleletype[1]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}	
		if(s1Alleletype[0]==s2Alleletype[1])
		{
			if(s1Alleletype[1]==s2Alleletype[0]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		output <<genotypes[sample1].positions[chrindex][sitei]<<"\t"<<numMatchMax<<endl;
	}
}

void comparison_usenamefilter_output(string sample1, string sample2, int chr, map<string, genotype> genotypes, ofstream &output, vector<string> filter)
{
	int chrindex=chr-1;
	for(int sitei=0;sitei<genotypes[sample1].positions[chrindex].size();sitei++)
	{
		if(find(filter.begin(),filter.end(),genotypes[sample1].sitenames[chrindex][sitei])==filter.end()) continue;
		int numMatchMax=0;
		string s1Alleletype=genotypes[sample1].alleletypes[chrindex][sitei];
		string s2Alleletype=genotypes[sample2].alleletypes[chrindex][sitei];


		if(s1Alleletype.find("-")!=string::npos || s2Alleletype.find("-")!=string::npos) continue;
		if(s1Alleletype.size()==0 || s1Alleletype.size()==0) continue;
		if(s1Alleletype[0]==s2Alleletype[0])
		{
			if(s1Alleletype[1]==s2Alleletype[1]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}	
		if(s1Alleletype[0]==s2Alleletype[1])
		{
			if(s1Alleletype[1]==s2Alleletype[0]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		output <<genotypes[sample1].positions[chrindex][sitei]<<"\t"<<numMatchMax<<endl;
	}
}

void comparison_usenamefilter_stats(string sample1, string sample2, int chr, map<string, genotype> genotypes, map<int, int> &matches, vector<string> filter)
{

	int chrindex=chr-1;
	//cout<<"1. "<<genotypes[sample1].alleletypes.size()<<"\t"<<genotypes[sample2].alleletypes.size()<<endl;
	//cout<<"2. "<<genotypes[sample1].alleletypes[chrindex].size()<<"\t"<<genotypes[sample2].alleletypes[chrindex].size()<<endl;
	for(int sitei=0;sitei<genotypes[sample1].positions[chrindex].size();sitei++)
	{

		if(find(filter.begin(),filter.end(),genotypes[sample1].sitenames[chrindex][sitei])==filter.end()) continue;

		int numMatchMax=0;
		string s1Alleletype=genotypes[sample1].alleletypes[chrindex][sitei];
		string s2Alleletype=genotypes[sample2].alleletypes[chrindex][sitei];


		if(s1Alleletype.find("-")!=string::npos || s2Alleletype.find("-")!=string::npos) continue;

		if(s1Alleletype[0]==s2Alleletype[0])
		{
			if(s1Alleletype[1]==s2Alleletype[1]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[1] && numMatchMax < 1) numMatchMax=1;
		}	
		if(s1Alleletype[0]==s2Alleletype[1])
		{
			if(s1Alleletype[1]==s2Alleletype[0]) numMatchMax=2;
			else if(s1Alleletype[1]!=s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		else
		{
			if(s1Alleletype[1]==s2Alleletype[0] && numMatchMax < 1) numMatchMax=1;
		}
		matches[genotypes[sample1].positions[chrindex][sitei]]=numMatchMax;

	}

}

//determine 
double contingencyTest(vector<int> presentBlock,vector<int> testBlock) //Cochran-Armitage test
{
	vector<int> N1; //presentBlock categories (wiki)
	N1.resize(3);
	N1[0]=0; N1[1]=0; N1[2]=0;
	for(auto pb : presentBlock)
	{
		N1[pb]++;
	}

	vector<int> N2; //testBlock categories
	N2.resize(3);
	N2[0]=0; N2[1]=0; N2[2]=0;
	for(auto tb : testBlock)
	{
		N2[tb]++;
	}

	int R1=N1[0]+N1[1]+N1[2];
	int R2=N2[0]+N2[1]+N2[2];
	//cout<<R1<<"\t"<<R2<<"\t"<<N1[0]<<"\t"<<N1[1]<<"\t"<<N1[2]<<"\t"<<N2[0]<<"\t"<<N2[1]<<"\t"<<N2[2]<<"\t";

	int C[3];
	C[0]=N1[0]+N2[0];
	C[1]=N1[1]+N2[1];
	C[2]=N1[2]+N2[2];
	int N=R1+R2;
	
	int t[3]={1,1,1};
	
	//statistic
	double T=0;
	for(int i=0; i<=2; i++)
	{
		T=T+t[i]*(R2*N1[i]-R1*N2[i]);
		//cout<<R2*N1[i]<<"\t"<<R1*N2[i]<<"\t";
	}
	

	double VarT=0;
	double VarT1=0;
	double VarT2=0;

	for(int i=0;i<=2;i++)
	{
		VarT1=VarT1+t[i]*t[i]*C[i]*(N-C[i]);
	}
	for(int i=0;i<=1;i++)
	{
		for(int j=i+1;j<=2;j++)
		{
			VarT2=VarT2+t[i]*t[j]*C[i]*C[j];
		}
	}

	VarT=R1*R2*(VarT1-2*VarT2)/T;
	
	double significanceTest=0;
	significanceTest=T/pow(VarT,0.5);

	//use negative value to get p
	double p=0.0;
	if(significanceTest>0) significanceTest=-significanceTest;
	p=phi(significanceTest);
	//cout<<T<<"\t"<<VarT<<"\t"<<significanceTest<<"\t"<<p<<endl;
	return p;
}

double gamma(double Z)
{
    const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
    const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI

    double D = 1.0 / (10.0 * Z);
    D = 1.0 / ((12 * Z) - D);
    D = (D + Z) * RECIP_E;
    D = pow(D, Z);
    D *= sqrt(TWOPI / Z);
 
    return D;
} 


static double igf(double S, double Z)
{
    if(Z < 0.0)
    {
	return 0.0;
    }
    double Sc = (1.0 / S);
    Sc *= pow(Z, S);
    Sc *= exp(-Z);
 
    double Sum = 1.0;
    double Nom = 1.0;
    double Denom = 1.0;
 
    for(int I = 0; I < 200; I++)
    {
	Nom *= Z;
	S++;
	Denom *= S;
	Sum += (Nom / Denom);
    }
 
    return Sum * Sc;
}

double chisqr(int Dof, double Cv)
{
    //cout<<"chisquare"<<endl;
    //cout<<Cv<<" "<<Dof<<endl;
    if(Cv < 0 || Dof < 1)
    {
        return 0.0;
    }
    double K = ((double)Dof) * 0.5;
    double X = Cv * 0.5;
    if(Dof == 2)
    {
	return exp(-1.0 * X);
    }

    double PValue = igf(K, X);
    //cout<<PValue<<endl;
    if(std::isnan(PValue) || std::isinf(PValue) || PValue <= 1e-8)
    {
        return 1e-14;
    } 

    PValue /= gamma(K);
    //cout<<PValue<<"\t"<<K<<"\t"<<gamma(K)<<endl;
    //PValue /= tgamma(K); 
	
    return (1.0 - PValue);
}


double contingencyTest1(vector<int> presentBlock,vector<int> testBlock) //chi-square test
{
	
	vector<double> N1; //presentBlock categories (wiki)
	N1.resize(3);
	N1[0]=0; N1[1]=0; N1[2]=0;
	for(auto pb : presentBlock)
	{
		N1[pb]++;
	}

	vector<double> N2; //testBlock categories
	N2.resize(3);
	N2[0]=0; N2[1]=0; N2[2]=0;
	for(auto tb : testBlock)
	{
		N2[tb]++;
	}

	double R1=N1[0]+N1[1]+N1[2];
	double R2=N2[0]+N2[1]+N2[2];

	double C[3];
	C[0]=N1[0]+N2[0];
	C[1]=N1[1]+N2[1];
	C[2]=N1[2]+N2[2];
	double N=R1+R2;

	double E1[3],E2[3];
	for(int i=0; i<3; i++)
	{
		E1[i]=C[i]*R1/(R1+R2);
		E2[i]=C[i]*R2/(R1+R2);
	}
	//cout<<"CV calc:"<<endl;
	//cout<<N1[0]<<" "<<N1[1]<<" "<<N1[2]<<" "<<N2[0]<<" "<<N2[1]<<" "<<N2[2]<<endl;
	//cout<<E1[0]<<" "<<E1[1]<<" "<<E1[2]<<" "<<E2[0]<<" "<<E2[1]<<" "<<E2[2]<<endl;	
	double chisquare=0;


	int df=2;	
	for(int i=0; i<3; i++)
	{
		if(E1[i]==0 || E2[i]==0 ) df--;
		else chisquare=chisquare+(N1[i]-E1[i])*(N1[i]-E1[i])/E1[i]+(N2[i]-E2[i])*(N2[i]-E2[i])/E2[i];
	}



	double p=chisqr(df, chisquare);
	return p;
}
	
	



void memberShip(map<int, int> matches, map<int, int> &membership, int windowSize, int stepSize, int minBlockSize)
{
	int presentBlockNumber=0; //any tested block are compared to the presentBlock, if not different, all sites in the block are assigned to the presentBlock
	vector<int> presentBlock;
	if(windowSize > matches.size()) cout<<"not enough sites"<<endl;
	else cout<<"enough sites"<<endl;

	for(auto mi=matches.begin();distance(matches.begin(),mi)<windowSize;mi++)
	{
		presentBlock.push_back(mi->second);
	}

	int isFirstDifferent=1;
	auto point1save=matches.begin();
	auto point2save=matches.begin();
	for(auto point1=matches.begin();;advance(point1,stepSize))
	{
		auto point2=point1;

		if(distance(point1,matches.end())<windowSize+windowSize) //deal with the last window situation
		{
			point2=matches.end();
			--point2;
		}
		else 
		{
			point2=point1; //point1 and point2-1 specify the block being tested
			advance(point2,windowSize);
		}
		//termination check

		//the test vector
		vector<int> testBlock;
		for(auto mi=point1;mi!=point2;mi++)
		{
			testBlock.push_back(mi->second);
		}

		//test if different
		double p=contingencyTest1(presentBlock,testBlock);

		//if different, update presentBlock
		//cout<<point1->first<<"\t"<<point2->first<<"\t"<<presentBlockNumber<<"\t"<<p<<endl;
		if(p<0.05 && p!=0)
		{
			if(isFirstDifferent==1)
			{
				isFirstDifferent=0;
				point1save=point1;
				point2save=point2;
			}
			else
			{

				isFirstDifferent=1;
				presentBlockNumber++;
				presentBlock=testBlock;
				//assign last block membership
				for(auto mi=point1save;mi!=point2save;mi++)
				{
					membership[mi->first]=presentBlockNumber;
				}
				//assign testblock membership
				for(auto mi=point1;mi!=point2;mi++)
				{
					membership[mi->first]=presentBlockNumber;
				}
			}
		}
		else
		{
			if(isFirstDifferent==0) 
			{
				//assign last block membership
				for(auto mi=point1save;mi!=point2save;mi++)
				{
					membership[mi->first]=presentBlockNumber;
				}
				isFirstDifferent=1;
			}

			//assign membership to the testblock			
			for(auto mi=point1;mi!=point2;mi++)
			{
				membership[mi->first]=presentBlockNumber;
			}
		}

		//cout<<"\t"<<point1->first<<"\t"<<point2->first<<"\t"<<presentBlockNumber<<"\t"<<p<<endl;
		

		if(distance(point2,matches.end())==1) break;
	}
	

	//If a block is too small (int terms of index), then probabily a transition state, assign -1.
	auto presentBlockStart=membership.begin();
	presentBlockNumber=0;

	/*
	for(auto mi=membership.begin();mi!=membership.end();mi++)
	{
		int testBlockNumber=mi->second;
		if(testBlockNumber!=presentBlockNumber)
		{
			//---see if block size < minBlockSize---
			if(distance(presentBlockStart,mi) < minBlockSize-1) 
			{
				for(auto tempi=presentBlockStart; tempi!=mi; tempi++)
				{
					tempi->second=-1;
				}
			}

			//---update presentBlockStart and presentBlockNumber
			presentBlockStart=mi;
			presentBlockNumber=mi->second;
		}
	}

	//check last block
	if(distance(presentBlockStart,membership.end())< minBlockSize-1)
	{
		for(auto tempi=presentBlockStart; tempi!=membership.end(); tempi++)
		{
			tempi->second=-1;
		}		
	}
	*/
}
	
class Block
{
public:
  int isTerminus;
  int length;
  int siteCounts;
  int chr;
  double meanMatches;
  int startBP;
};

void extractBlockInformation(map<int,int> matches, map<int,int> membership, vector<Block> &blocks)
{
  int presentblockIndex=-1;
  int presentBlockStartBP=0;
  for(auto record : membership)
  {
      if(record.second!=presentblockIndex)
      {
	presentblockIndex=record.second;
	Block emptyblock;
	blocks.push_back(emptyblock);
	blocks.back().siteCounts=0;
	presentBlockStartBP=record.first;
	blocks.back().startBP=record.first;
      }
      blocks.back().siteCounts++; 
      blocks.back().length=record.first-presentBlockStartBP;
      blocks.back().meanMatches+=matches[record.first];
  }

  for(auto bi=blocks.begin();bi!=blocks.end();bi++)
  {
    (bi->meanMatches)=(bi->meanMatches)/(bi->siteCounts);
  }
}

string combineString(string s1, string s2)
{
	string s;
	if(s1>s2)
	{
		s=s2+"\t"+s1;
	}
	else
	{
		s=s1+"\t"+s2;
	}

return s;
}

void forceInsertRecord(string s1, string s2, double d, map<string, double> &m) //insert relationship d of individuals s1 and s2, into m
{
	if((s1!="") && (s2!=""))
	{
		string s=combineString(s1,s2);
		m[s]=d;
	}
}

int main(int argc, char * argv[])
{
	string pairlist=argv[5];

	string variablelistfile=argv[2];
	string genotypefile=argv[1];
	string snpstatefile=argv[3];
	string samplelistfile=argv[4];
	string pairout=argv[6];
	string isskip="0";
// 	
	string individual1;
	string individual2;

	map<string, location> snpLocations;
	map<string, genotype> genotypes;

	ifstream input1(variablelistfile.c_str());
	ifstream input2(genotypefile.c_str());
	ifstream input3(snpstatefile.c_str());

	cout<<"---Parsing bp to cm mapping..."<<endl;
	ParseMap(input3, snpLocations);
	cout<<"---DONE"<<endl;
	cout<<endl;

	cout<<"---Reading in SNP profiles..."<<endl;
	ParseInput(input2,  snpLocations, genotypes);
	cout<<"---DONE"<<endl;
	cout<<endl;
	
	//read in the sites that have MAF >=0.4 in unrelated subjects
	cout<<"---Reading in high MAF site list..."<<endl;
	vector<string> filter;
	string line;
	for(;getline(input1,line);) filter.push_back(line);
	cout<<"---DONE"<<endl;
	cout<<endl;

	//a list of samples
	cout<<"---Reading in sample list..."<<endl;
	ifstream sampleList(samplelistfile.c_str());
	vector<string> samples;
	for(;getline(sampleList,line);)
	{
	  samples.push_back(line);
	}
	cout<<"---DONE"<<endl;
	cout<<endl;



	//a list of paired samples
	cout<<"---Reading in sample relations..."<<endl;
	ifstream paired(pairlist.c_str());
	map<string,double> sampleRelation;
	for(;getline(paired,line);)
	{
		istringstream linestream(line);
		string f1,f2,f3;
		getline(linestream,f1,'\t');
		getline(linestream,f2,'\t');		
		getline(linestream,f3,'\t');
		forceInsertRecord(f1, f2, stod(f3), sampleRelation);
	}
	cout<<"---DONE"<<endl;
	cout<<endl;


	//---centromere position---
	vector<double> centromere(22);
	centromere[0]=125000000;
	centromere[1]=93300000;
	centromere[2]=91000000;
	centromere[3]=50400000;
	centromere[4]=48400000;
	centromere[5]=61000000;
	centromere[6]=59900000;
	centromere[7]=45600000;
	centromere[8]=49000000;
	centromere[9]=40200000;
	centromere[10]=53700000;
	centromere[11]=35800000;
	centromere[12]=17900000;
	centromere[13]=17600000;
	centromere[14]=19000000;
	centromere[15]=36600000;
	centromere[16]=24000000;
	centromere[17]=17200000;
	centromere[18]=26500000;
	centromere[19]=27500000;
	centromere[20]=13200000;
	centromere[21]=14700000;
	

	//compare between all pairs of samples
	cout<<"---Starting analysis core..."<<endl;
	ofstream blockstats(pairout.c_str());


	for(int sample1i=0;sample1i<samples.size();sample1i++)
	{
	    for(int sample2i=sample1i+1;sample2i<samples.size();sample2i++)
	    {

		
		individual1=samples[sample1i];
		individual2=samples[sample2i];

		//check if relationship between the two people are established
		
		string s=combineString(individual1, individual2);
		double relationship;
		if(sampleRelation.find(s)==sampleRelation.end()) continue;
		else relationship=sampleRelation[s];

		cout<<"------Testing "<<individual1<<" and "<<individual2<<"..."<<endl;
		
		vector<Block> allChromosomeBlocks;

		for(int chri=1; chri<=22; chri++)
		{
		    //-----extract number of matches----------
		    if(genotypes[individual1].positions[chri-1].size()==0 || genotypes[individual1].positions[chri-1].size()==0) continue;
		    map<int, int> matches; // for a given chr and a pair of samples, position (bp) -> # of max matches
		    comparison_usenamefilter_stats(individual1, individual2, chri, genotypes, matches, filter);

		    //-----compute block membership----------
		    map<int, int> membership;
		    cout<<chri<<endl;
		    memberShip(matches, membership, 50, 50, 100);
		    //cout<<matches.size()<<" "<<membership.size()<<endl;
		    //-----extract block information---------
		    vector<Block> blocks;
		    extractBlockInformation(matches,membership,blocks);
		    for(int bli=0;bli<blocks.size();bli++) blocks[bli].chr=chri;
			allChromosomeBlocks.insert(allChromosomeBlocks.end(), blocks.begin(), blocks.end());
		}

		//-----output------
		//1. output all  matches
		for(int bi=0;bi<allChromosomeBlocks.size();bi++)
		{
			if(isskip=="1") //skip blocks with centromere
			{
				if((allChromosomeBlocks[bi].startBP<=centromere[allChromosomeBlocks[bi].chr-1]) && (allChromosomeBlocks[bi].startBP+allChromosomeBlocks[bi].length>=centromere[allChromosomeBlocks[bi].chr-1])) continue;
			}
			blockstats<<individual1<<" "<<individual2<<"\t"<<allChromosomeBlocks[bi].length<<"\t"<<allChromosomeBlocks[bi].chr<<"\t"<<allChromosomeBlocks[bi].meanMatches<<"\t"<<relationship<<"\t"<<allChromosomeBlocks[bi].startBP<<endl;
		}

		cout<<"------DONE"<<endl;
	    }
	}
	cout<<"---ALL DONE"<<endl;
}
