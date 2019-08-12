// g++ -std=c++11 deserialize.cpp -o deserialize.o -I /igm/projects/kmer_dragen/parallel-hashmap/parallel_hashmap -I /igm/projects/kmer_dragen/cereal/include -O3

//  Example execution
// ---------------
// kmer_ref_segmentak.cereal
// kmer_ref_segmental.cereal
// qsub -pe smp 14 -b y -j y -cwd /igm/projects/kmer_dragen/exp_code/ssw_multi.o /igm/projects/kmer_dragen/exp_code/kmer_ref_segmentak.cereal /igm/projects/kmer_dragen/exp_code/kmer_ref_segmental.cereal /igm/projects/kmer_dragen/akj_father_HG003/test_reads.txt /igm/projects/kmer_dragen/akj_father_HG003/test_reads.txt /igm/projects/kmer_dragen/akj_father_HG003/test_reads.txt /igm/projects/kmer_dragen/exp_code/output/parallel_0.txt /igm/projects/kmer_dragen/exp_code/output/parallel_1.txt /igm/projects/kmer_dragen/exp_code/output/parallel_2.txt


#include <array>
#include <queue>
#include <thread>
#include <iostream>
#include <bitset>
#include <cinttypes>
#include "parallel_hashmap/phmap.h"
#include "cereal/types/unordered_map.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/bitset.hpp"
#include "cereal/archives/binary.hpp"
#include <fstream>
#include <random>
#include <chrono>
#include <functional>
#include <cstdio>
#include <future>
#include <iostream>
#include <vector>
#include "ssw.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <queue>

using phmap::flat_hash_map;
using namespace std;
int ksize = 21;


// char complement(char n)
// {   
//     switch(n)
//     {   
//     case 'A':
//         return 'T';
//     case 'T':
//         return 'A';
//     case 'G':
//         return 'C';
//     case 'C':
//         return 'G';
//     case 'N':
//         return 'N';        
//     }   
//     assert(false);
//     return ' ';
// }  

std::vector<int> topK(int indices_tmp[5], int k)
{
  std::priority_queue< std::pair<double, int>, std::vector< std::pair<double, int> >, std::greater <std::pair<double, int> > > q;
  for (int i = 0; i < 5; ++i) 
  {
    if(q.size()<k)
        q.push(std::pair<double, int>(indices_tmp[i], i));
    else if(q.top().first < indices_tmp[i])
    {
        q.pop();
        q.push(std::pair<double, int>(indices_tmp[i], i));
    }
  }
  k = q.size();
  std::vector<int> res(k);
  for (int i = 0; i < k; ++i) {
    res[k - i - 1] = q.top().second;
    q.pop();
  }
  return res;
}


int ReadMapAlign(string &fastq, string &fileout, const flat_hash_map<bitset<42>, long>& rev_strand_kmers, const flat_hash_map<bitset<42>, long>& fwd_strand_kmers, string &ref_str_f, string &ref_str_r)
{
    
    int32_t l, m, k, match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;
    
    static const int8_t nt_table[128] = 
    {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };    

    string output;
    // int reverse_counts = 0;
    // int reverse_counts_abv_z = 0;
    // int forward_counts = 0;
    // int forward_counts_abv_z = 0;
    // load and parse fastq
    // -----------    
    // collect the hashmap and then run the analysis and output what we need
    std::ofstream fout(fileout);
    ifstream input_stream(fastq);
    string read;
    int c = 0; //track the location of the read with respect to the file
    while (getline(input_stream, read)) 
    {
        if (c == 4) 
        {
            c = 0;
        }

        if (c == 1)
        {
            std::transform(read.begin(), read.end(), read.begin(), ::toupper);
            int read_len = read.length();
            string kmer;
            int max_i[2] = {0,0}; // store the max value from strand 1 and 0
            // completement
            // -------            
            int g = 0; // keep track of how many pos. we placed
            long ref_pos[5] = {0,0,0,0,0}; // only store max 5 values
            int kmer_pos[5] = {0,0,0,0,0}; // only store max 5 values
            int max_count[5] = {0,0,0,0,0}; // only store max 5 values
            // reverse
            // -------
            int g_r = 0; // keep track of how many pos. we placed
            long ref_pos_r[5] = {0,0,0,0,0}; // only store max 5 values
            int kmer_pos_r[5] = {0,0,0,0,0}; // only store max 5 values
            int max_count_r[5] = {0,0,0,0,0}; // only store max 5 values                        
            // int r = 0;            
            // for (auto& table : fmap)
            for (int r = 0; r <2; ++r)
            {

                // flat_hash_map<bitset<42>, long> table = f.get();
                //  sliding window
                for (int i = 0; i < read_len-ksize-1; i+=2) // increment by 2
                {
                    bitset<42> bs_kmer;
                    int k_p = 0; // bits need to be counted by 2s
                    // convert to bitset
                    try
                    {
                        kmer = read.substr(i,ksize);
                    }
                    catch (const std::out_of_range& e)
                    {
                        cout << "Out of Range error.";
                    }
                    if ((kmer.find('N') == string::npos) && !kmer.empty())
                    {
                        for (int j = 0; j < ksize; ++j)
                        {
                            if (kmer[j] == 'C')
                            {
                                bs_kmer.set(k_p+1);
                            }
                            if (kmer[j] == 'T')
                            {
                                bs_kmer.set(k_p+1);
                                bs_kmer.set(k_p);
                            }
                            if (kmer[j] == 'G')
                            {
                                bs_kmer.set(k_p);
                            } 
                            k_p += 2;
                        }
                        // extract the reference position where a kmer_matches
                        if ((fwd_strand_kmers.find(bs_kmer) != fwd_strand_kmers.end()) && (r == 0))
                        {
                            auto kmer_it = fwd_strand_kmers.find(bs_kmer);
                            long n_ref_pos = kmer_it->second;
                            if (g == 5)
                            {
                                break;
                            }                        
                            if ((ref_pos[0] == 0) && (g == 0))
                            {
                                ref_pos[0] = n_ref_pos;
                                kmer_pos[0] = i;
                                max_count[0] += 1;
                            }
                            else
                            {
                                for (int w=0; w<g+1; ++w)
                                {
                                    if ((abs(n_ref_pos - ref_pos[w]) <= read_len))
                                    {
                                        max_count[w] += 1;
                                    }
                                    else if (((n_ref_pos - ref_pos[w]) > read_len) && (w < 5)) 
                                    {
                                        ref_pos[w+1] = n_ref_pos;
                                        kmer_pos[w+1] = i;
                                        g += 1;
                                    }
                                }
                            }
                        }
                        else if ((rev_strand_kmers.find(bs_kmer) != rev_strand_kmers.end()) && (r == 1))
                        {
                            auto kmer_it = rev_strand_kmers.find(bs_kmer);
                            long n_ref_pos = kmer_it->second;
                            if (g_r == 5)
                            {
                                break;
                            }
                            if ((ref_pos_r[0] == 0) && (g_r == 0))
                            {
                                ref_pos_r[0] = n_ref_pos;
                                kmer_pos_r[0] = i;
                                max_count_r[0] += 1;
                            }
                            else
                            {
                                for (int w=0; w<g_r+1; ++w)
                                {
                                    if ((abs(n_ref_pos - ref_pos_r[w]) <= read_len))
                                    {
                                        max_count_r[w] += 1;
                                    }
                                    else if (((n_ref_pos - ref_pos_r[w]) > read_len) && (w < 5)) 
                                    {
                                        ref_pos_r[w+1] = n_ref_pos;
                                        kmer_pos_r[w+1] = i;
                                        g_r += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                // after computing the mapping metrics for the read assign a score to differentiate the mapping of the strand
                if (r == 0)
                {
                    // int map_max = 0;
                    for (int z=0; z<5; ++z)
                    {
                        if (max_count[z] > max_i[0])
                        {
                           // map_max = max_count[z]; // assign the number of mapped regions to the placeholder
                            max_i[0] = max_count[z];
                        }
                    }
                    // max_i[0] = map_max; // assign strand 0 the highest mapping count to max_i array.
                }
                else if (r == 1)
                {
                    // int map_max = 0;
                    for (int z=0; z<5; ++z)
                    {
                        if (max_count_r[z] > max_i[1])
                        {
                           // map_max = max_count_r[z];
                            max_i[1] = max_count_r[z];
                        }
                    }
                    // max_i[1] = map_max; // assign strand 1 the highest mapping count to max_i array.
                }
                // r += 1;
            }
            // cout << ": F "+ to_string(max_count[0]) +':'+to_string(max_count[1]) +':'+to_string(max_count[2])+':'+to_string(max_count[3])+':'+to_string(max_count[4])+ ": R "+ to_string(max_count_r[0]) +':'+to_string(max_count_r[1]) +':'+to_string(max_count_r[2])+':'+to_string(max_count_r[3])+':'+to_string(max_count_r[4])+'\n';
            // iterate over the reference positions and extract the reference sequence of interest
            // string ref_output = "abc";
            // string ref_output;
            // Forward strand check
            // ---------------------                        
            if (max_i[0] >= max_i[1]) // arbitrary bias towards strand 0
            {
                // forward_counts += 1;
                // int top_indices[2] = {0,0};
                std::vector<int>  top_indices = topK(max_count, 2);
                // for (int z=0; z<5; ++z)
                // {
                //     if (max_count[z] > top_indices[0])
                //     {
                //         top_indices[0] = z;
                //     }
                //     if ((max_count[z] < top_indices[0]) && (max_count[z] > top_indices[1]))
                //     {
                //         top_indices[1] = z;
                //     }
                // }
                string tmp_ref_output;
                string sw_scores;
                // string tmp_pos;
                for (int idx =0; idx <2; ++idx)
                {                    
                    // check for reverse max against the second highest forward max
                    if ((idx == 1) && (max_count[top_indices[1]] < max_i[1]))
                    {                        
                        // int index_r = 0;
                        // for (int z=0; z<5; ++z)
                        // {
                        //     if (max_count_r[z] > index_r)
                        //     {
                        //         index_r = z;
                        //     }
                        // }
                        int index_r = topK(max_count_r, 1)[0];
                        int val = max_count_r[index_r];
                        long p = ref_pos_r[index_r];
                        if (val > 0)
                        {
                            // reverse_counts_abv_z += 1;
                            long rp = p - kmer_pos_r[index_r];
                            try 
                            {
                                string ref_output = ref_str_r.substr(rp, read_len);

                                // SSW  STEP
                                // ------------------------------
                                s_profile* profile;
                                int8_t* num = (int8_t*)malloc(read_len);
                                int8_t* ref_num = (int8_t*)malloc(read_len);
                                s_align* result;
                                int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
                                for (l = k = 0; l < 4; ++l)
                                {
                                    for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;
                                        mat[k++] = 0;
                                }
                                for (m = 0; m < 5; ++m) mat[k++] = 0;
                                for (m = 0; m < read_len; ++m) num[m] = nt_table[(int)read[m]];    
                                profile = ssw_init(num, read_len, mat, 5, 2);
                                for (m = 0; m < read_len; ++m) ref_num[m] = nt_table[(int)ref_output[m]];
                                result = ssw_align(profile, ref_num, read_len, gap_open, gap_extension, 1, 0, 0, 15);
                                // uint16_t score = result->score1;
                                // sw_scores = "," + to_string(score);
                                // free(result->cigar);
                                // free(result);
                                profile = NULL;
                                result = NULL;
                                // align_destroy(result);
                                // init_destroy(profile);
                                free(mat);
                                free(ref_num);
                                free(num);                        
                            }
                            catch (const std::out_of_range& e) 
                            {
                               cout << "Out of Range error.";
                            }                
                        // tmp_ref_output += " @ " + to_string(idx) +": F " + to_string(max_count[top_indices[1]]) + ": R "+ to_string(max_count_r[0]) +':'+to_string(max_count_r[1]) +':'+to_string(max_count_r[2])+':'+to_string(max_count_r[3])+':'+to_string(max_count_r[4]) + " pos " + to_string(rp);
                            // tmp_ref_output += "," + to_string(idx)+":"+sw_scores;
                        }
                    }
                    else
                    {
                        int index_f = top_indices[idx];
                        long p = ref_pos[index_f];
                        int val = max_count[index_f];
                        // tmp_pos = ":"+to_string(max_count[top_indices[0]])+":"+to_string(max_count_r[top_indices[1]])+"@"+to_string(idx);
                        if (val > 0)
                        {
                            // forward_counts_abv_z += 1;
                            long rp = p - kmer_pos[index_f];
                            try 
                            {
                                string ref_output = ref_str_f.substr(rp, read_len);

                                // SSW  STEP
                                // ------------------------------
                                // s_profile* profile;
                                int8_t* num = (int8_t*)malloc(read_len);
                                int8_t* ref_num = (int8_t*)malloc(read_len);
                                // s_align* result;
                                int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
                                // for (l = k = 0; l < 4; ++l)
                                // {
                                    // for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;
                                        // mat[k++] = 0;
                                // }
                                // for (m = 0; m < 5; ++m) mat[k++] = 0;
                                // for (m = 0; m < read_len; ++m) num[m] = nt_table[(int)read[m]];    
                                // profile = ssw_init(num, read_len, mat, 5, 2);
                                // for (m = 0; m < read_len; ++m) ref_num[m] = nt_table[(int)ref_output[m]];
                                // result = ssw_align(profile, ref_num, read_len, gap_open, gap_extension, 1, 0, 0, 15);
                                // uint16_t score = result->score1;
                                // sw_scores = "," + to_string(score);
                                // align_destroy(result);
                                // init_destroy(profile);
                                // free(result);
                                // free(result->cigar);
                                free(mat);
                                free(ref_num);
                                free(num);

                            }
                            catch (const std::out_of_range& e) 
                            {
                               cout << "Out of Range error.";
                            }
                           // tmp_ref_output += " @ " + to_string(idx) +": F " + to_string(max_count[top_indices[idx]]) + ": F "+ to_string(max_count[0]) +':'+to_string(max_count[1]) +':'+to_string(max_count[2])+':'+to_string(max_count[3])+':'+to_string(max_count[4]) + ":pos " + to_string(rp);
                            // tmp_ref_output += "," + to_string(idx)+":"+sw_scores;
                        }
                    }
                    
                }
                if (!tmp_ref_output.empty())
                {
                    // output += (read + "," + tmp_ref_output +'\n');
                    // output += read + '\n';
                    // output += to_string(max_count[max_i]) + '\n';
                }                
            }
            // Reverse strand check
            // ---------------------
            else
            {
                // reverse_counts += 1;
                // int top_indices[2] = {0,0};
                std::vector<int>  top_indices = topK(max_count_r,2);
                // for (int z=0; z<5; ++z)
                // {
                //     if (max_count_r[z] > top_indices[0])
                //     {
                //         top_indices[0] = z;
                //     }
                //     if ((max_count_r[z] < top_indices[0]) && (max_count_r[z] > top_indices[1]))
                //     {
                //         top_indices[1] = z;
                //     }
                // }
                string tmp_ref_output;
                string sw_scores;
                for (int idx =0; idx <2; ++idx)
                {
                    if ((idx == 1) && (max_count_r[top_indices[1]] < max_i[0]))
                    {
                        // int f_index = 0;
                        // for (int z=0; z<5; ++z)
                        // {
                        //     if (max_count[z] > f_index)
                        //     {
                        //         f_index = z;
                        //     }
                        // }
                        int f_index = topK(max_count, 1)[0];
                        long p = ref_pos[f_index];
                        int val = max_count[f_index];
                        if (val > 0)
                        {
                            // forward_counts_abv_z += 1;
                            long rp = p - kmer_pos[f_index];
                            try 
                            {
                                string ref_output = ref_str_f.substr(rp, read_len);
                                // SSW  STEP
                                // ------------------------------
                                // s_profile* profile;
                                int8_t* num = (int8_t*)malloc(read_len);
                                int8_t* ref_num = (int8_t*)malloc(read_len);
                                // s_align* result;
                                int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
                                // for (l = k = 0; l < 4; ++l)
                                // {
                                //     for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;
                                //         mat[k++] = 0;
                                // }
                                // for (m = 0; m < 5; ++m) mat[k++] = 0;
                                // for (m = 0; m < read_len; ++m) num[m] = nt_table[(int)read[m]];    
                                // profile = ssw_init(num, read_len, mat, 5, 2);
                                // for (m = 0; m < read_len; ++m) ref_num[m] = nt_table[(int)ref_output[m]];
                                // result = ssw_align(profile, ref_num, read_len, gap_open, gap_extension, 1, 0, 0, 15);
                                // uint16_t score = result->score1;
                                // sw_scores = "," + to_string(score);
                                // align_destroy(result);
                                // init_destroy(profile);
                                // free(result);
                                // free(result->cigar);
                                free(mat);
                                free(ref_num);
                                free(num);

                            }
                            catch (const std::out_of_range& e) 
                            {
                               cout << "Out of Range error.";
                            }
                        // tmp_ref_output += " @ " + to_string(idx) +": F " + to_string(max_count[top_indices[idx]]) + ": F "+ to_string(max_count[0]) +':'+to_string(max_count[1]) +':'+to_string(max_count[2])+':'+to_string(max_count[3])+':'+to_string(max_count[4]) + ":pos " + to_string(rp);    
                            // tmp_ref_output += "," + to_string(idx)+":"+sw_scores;
                        }

                    }
                    // under usual circumstances
                    else
                    {
                        int r_index = top_indices[idx];
                        long p = ref_pos_r[r_index];
                        int val = max_count[r_index];
                        if (val > 0)
                        {
                            // reverse_counts_abv_z += 1;
                            long rp = p - kmer_pos_r[r_index];
                            try 
                            {
                                string ref_output = ref_str_r.substr(rp, read_len);

                                // SSW  STEP
                                // ------------------------------
                                // s_profile* profile;
                                int8_t* num = (int8_t*)malloc(read_len);
                                int8_t* ref_num = (int8_t*)malloc(read_len);
                                // s_align* result;
                                int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
                                // for (l = k = 0; l < 4; ++l)
                                // {
                                //     for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;
                                //         mat[k++] = 0;
                                // }
                                // for (m = 0; m < 5; ++m) mat[k++] = 0;
                                // for (m = 0; m < read_len; ++m) num[m] = nt_table[(int)read[m]];    
                                // profile = ssw_init(num, read_len, mat, 5, 2);
                                // for (m = 0; m < read_len; ++m) ref_num[m] = nt_table[(int)ref_output[m]];
                                // result = ssw_align(profile, ref_num, read_len, gap_open, gap_extension, 1, 0, 0, 15);
                                // uint16_t score = result->score1;
                                // sw_scores = "," + to_string(score);
                                // align_destroy(result);
                                // init_destroy(profile);
                                // free(result);
                                // free(result->cigar);
                                free(mat);
                                free(ref_num);
                                free(num);                        
                            }
                            catch (const std::out_of_range& e) 
                            {
                               cout << "Out of Range error.";
                            }
                        // tmp_ref_output += " @ " + to_string(idx) +": F " + to_string(max_count[top_indices[1]]) + ": R "+ to_string(max_count_r[0]) +':'+to_string(max_count_r[1]) +':'+to_string(max_count_r[2])+':'+to_string(max_count_r[3])+':'+to_string(max_count_r[4]) + " pos " + to_string(rp);                   
                            // tmp_ref_output += "," + to_string(idx)+":"+sw_scores;
                        }    
                    }
                    
                }
                if (!tmp_ref_output.empty())
                {
                    // output += (read + "," + tmp_ref_output +'\n');
                    // output += read + '\n';
                    // output += to_string(max_count[max_i]) + '\n';
                }                
            }
        }
        c += 1;
    }
    fout << output;
    // fout << c << ' ' << forward_counts << ' ' << forward_counts_abv_z << ' ' << reverse_counts << ' ' << reverse_counts_abv_z;
    return 0;
}



// build the dictionaries
// -----------------
flat_hash_map<bitset<42>, long> getTable(string &fp)
{
    using MapType = phmap::flat_hash_map<bitset<42>, long>;
    MapType table;
    ifstream is(fp, ios::binary);
    cereal::BinaryInputArchive archive_in(is);
    size_t table_size;
    archive_in(table_size);
    table.reserve(table_size);
    archive_in(table);
    return table;
}

string getFullSequence(string &fp)
{
    std::ifstream reference_file(fp);
    std::string ref_str;
    reference_file.seekg(0, std::ios::end);
    ref_str.reserve(reference_file.tellg());
    reference_file.seekg(0, std::ios::beg);
    ref_str.assign((std::istreambuf_iterator<char>(reference_file)),
            std::istreambuf_iterator<char>());
    return ref_str;

}

int main(int argc, char* argv[])
{   

    string full_seq_fp[2] = {"path_to_forward", "path_to_reverse"};
    
    vector<future<string>> futures_full_seq;

    for(int i = 0; i < 2; ++i)
    {
        futures_full_seq.emplace_back(std::async(std::launch::async, getFullSequence, std::ref(full_seq_fp[i])));
    }
    string ref_str_fwd = futures_full_seq[0].get();
    string ref_str_rev = futures_full_seq[1].get();
    // long ref_len = ref_str_rev.length();
    // Load the whole reference into memory
    // ------------------------------------
    // std::ifstream reference_file("/igm/home/axr102/smith_waterman/human_g1k_v37_decoy.str");
    // std::string ref_str;
    // reference_file.seekg(0, std::ios::end);
    // ref_str.reserve(reference_file.tellg());
    // reference_file.seekg(0, std::ios::beg);
    // ref_str.assign((std::istreambuf_iterator<char>(reference_file)),
    //         std::istreambuf_iterator<char>());
    // ref_str.erase(std::remove(ref_str.begin(), ref_str.end(), '\n'), ref_str.end());
    // long ref_len = ref_str.length();

    // vector of files
    string path_to_table[2] = {argv[1], argv[2]};
    // deserialize tables
    // -----------
    vector<future<flat_hash_map<bitset<42>, long>>> futures;

    for(int i = 0; i < 2; ++i)
    {
        futures.emplace_back(std::async(std::launch::async, getTable, std::ref(path_to_table[i])));
    }
    
    // get dictionaries 
    const flat_hash_map<bitset<42>, long> fwd_strand_kmers = futures[0].get();
    const flat_hash_map<bitset<42>, long> rev_strand_kmers = futures[1].get();
    futures.clear();

    string fastq[3] = {argv[3], argv[4],argv[5]};
    string outputs[3] = {argv[6], argv[7],argv[8]};

    vector<future<int>> next_step;

    // placeholder for file division
    for(int i = 0; i < 3; ++i)
    {
        next_step.emplace_back(std::async(std::launch::async, ReadMapAlign, std::ref(fastq[i]), std::ref(outputs[i]), std::ref(rev_strand_kmers), std::ref(fwd_strand_kmers), std::ref(ref_str_fwd), std::ref(ref_str_rev)));
    }

    return 0;
}
