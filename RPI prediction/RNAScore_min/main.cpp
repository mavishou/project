#include "iostream"
#include "fstream"
#include "RNA.h"
#include "protein.h"

using namespace std;

int main( int argc , char *argv[] )
{
    string filename = argv[2];
    protein *testp = new protein ;
    testp->load_file( filename );
    testp->gen_score1();
    testp->gen_score2();
    testp->gen_score3();
    testp->gen_score4();
    testp->gen_score5();

    filename = argv[1];
    RNA *testr = new RNA;
    testr->init( filename );
//  cout << "step 1 is running..." << endl;
    testr->stepFirst();
//    cout << "step 2 is running..." << endl;
    testr->stepSecond();
//    cout << "step 3 is running..." << endl;

    testr->humanprotein_1 = testr->load_humanprotein( testp->score1 , testr->humanprotein_1 , testr->len_humanprotein_1 );
    testr->humanprotein_2 = testr->load_humanprotein( testp->score2 , testr->humanprotein_2 , testr->len_humanprotein_2 );
    testr->humanprotein_3 = testr->load_humanprotein( testp->score3 , testr->humanprotein_3 , testr->len_humanprotein_3 );
    testr->humanprotein_4 = testr->load_humanprotein( testp->score4 , testr->humanprotein_4 , testr->len_humanprotein_4 );
    testr->humanprotein_5 = testr->load_humanprotein( testp->score5 , testr->humanprotein_5 , testr->len_humanprotein_5 );

    testr->stepThird();
    testr->protein_name = testp->all_name;
    //cout << "output as file \"score.txt\"..." << endl;
    testr->output();

    return 0;
}
