#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <functional>
#include <random>
#include <bitset>
#include <cassert>
#include <iomanip>
#include <algorithm>



uint32_t mask_forward = 0x80c18000;
uint32_t mask_backward = 0x160a0e06;
uint32_t mask_match = 0x07000003;
int d_f = 5, d_b = 10, d_m = 5;

static constexpr uint32_t T(int j) {
    return (j >= 0 && j <= 15) ? 0x79cc4519U : 0x7a879d8aU;
}

uint32_t rol(uint32_t x, int n) {
    n &= 31;
    if (n == 0) return x;
    return (x << n) | (x >> (32 - n));
}

uint32_t ror(uint32_t x, int n) {
    n &= 31;
    if (n == 0) return x;
    return (x >> n) | (x << (32 - n));
}

uint32_t XOR(uint32_t x, uint32_t y) { return x ^ y; }

uint32_t P0(uint32_t X) { return X ^ rol(X, 9) ^ rol(X, 17); }
uint32_t P1(uint32_t X) { return X ^ rol(X, 15) ^ rol(X, 23); }

uint32_t InvP0(uint32_t X) {
    return X ^ rol(X, 2) ^ rol(X, 3) ^ rol(X, 9) ^ rol(X, 11) ^ rol(X, 17) ^ rol(X, 18) ^ rol(X, 19) ^ rol(X, 27);
}
uint32_t InvP1(uint32_t X) {
    return X ^ rol(X, 5) ^ rol(X, 13) ^ rol(X, 14) ^ rol(X, 15) ^ rol(X, 21) ^ rol(X, 23) ^ rol(X, 29) ^ rol(X, 30);
}

uint32_t FF(int j, uint32_t X, uint32_t Y, uint32_t Z) {
    if (j >= 0 && j <= 15) return X ^ Y ^ Z;
    return (X & Y) | (X & Z) | (Y & Z);
}

uint32_t GG(int j, uint32_t X, uint32_t Y, uint32_t Z) {
    if (j >= 0 && j <= 15) return X ^ Y ^ Z;
    return (X & Y) | ((0xFFFFFFFFU ^ X) & Z);
}




void generateRandomUInt32Array(uint32_t* array, size_t size) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<uint32_t> dis(0, UINT32_MAX); // 0 到 2^32 - 1

    for (size_t i = 0; i < size; ++i) {
        array[i] = dis(gen);
    }
}

void StepFunction(uint32_t state[8], uint32_t w, uint32_t w_prime, uint8_t iR) {
    uint32_t ff, gg, ss1, ss2, tt1, tt2, newB, newF, newE;

    ff = FF(iR, state[0], state[1], state[2]);      //多数函数
    gg = GG(iR, state[4], state[5], state[6]);      //选择函数

    ss1 = rol(rol(state[0], 12) + state[4] + rol(T(iR), (iR % 32)), 7);
    ss2 = XOR(ss1, rol(state[0], 12));

    tt1 = ff + state[3] + ss2 + w_prime;
    tt2 = gg + state[7] + ss1 + w;

    newB = rol(state[1], 9);
    newF = rol(state[5], 19);
    newE = P0(tt2);

    state[7] = state[6];            //H <- G
    state[6] = newF;                //G <- F <<< 19
    state[5] = state[4];            //F <- E
    state[4] = newE;                //E <- P0
    state[3] = state[2];            //D <- C
    state[2] = newB;                //C <- B <<< 9
    state[1] = state[0];            //B <- A
    state[0] = tt1;                 //A <- tt1 
}

void StepFunctionSP(uint32_t state[8], uint8_t iR) {
    uint32_t ff, gg, ss1, ss2, tt1, tt2, newB, newF, newE;

    ff = FF(iR, state[0], state[1], state[2]);      //多数函数
    gg = GG(iR, state[4], state[5], state[6]);      //选择函数

    ss1 = rol(rol(state[0], 12) + state[4] + rol(T(iR), (iR % 32)), 7);
    ss2 = XOR(ss1, rol(state[0], 12));

    tt1 = ff + state[3] + ss2;
    tt2 = gg + state[7] + ss1;

    newB = rol(state[1], 9);
    newF = rol(state[5], 19);
    newE = P0(tt2);

    state[7] = state[6];            //H <- G
    state[6] = newF;                //G <- F <<< 19
    state[5] = state[4];            //F <- E
    state[4] = newE;                //E <- P0
    state[3] = state[2];            //D <- C
    state[2] = newB;                //C <- B <<< 9
    state[1] = state[0];            //B <- A
    state[0] = tt1;                 //A <- tt1 
}

void InvStepFunction(uint32_t state[8], uint32_t w, uint32_t w_prime, uint8_t iR) {
    uint32_t ff, gg, ss1, ss2, tt1, tt2, newD, newH, newE;

    ff = FF(iR, state[1], ror(state[2], 9), state[3]);
    ss1 = rol(rol(state[1], 12) + state[5] + rol(T(iR), (iR % 32)), 7);
    ss2 = XOR(ss1, rol(state[1], 12));

    gg = GG(iR, state[5], ror(state[6], 19), state[7]);

    newD = state[0] - w_prime - ff - ss2;
    newH = InvP0(state[4]) - w - gg - ss1;

    state[0] = state[1];            // A <- B
    state[1] = ror(state[2], 9);            // B <- C >>> 9
    state[2] = state[3];            // C <- D
    state[3] = newD;                // D <- newD
    state[4] = state[5];            // E <- F
    state[5] = ror(state[6], 19);            // F <- G >>> 19
    state[6] = state[7];            // G <- H
    state[7] = newH;                // H <-newH
}





template <typename KeyType, typename ValueType>
class MultiHashTable {
private:
    struct HashNode {
        KeyType key;
        ValueType value;
        HashNode(const KeyType& k, const ValueType& v) : key(k), value(v) {}
    };

    std::vector<std::list<HashNode>> table;
    size_t size;
    size_t capacity;

    size_t hash(const KeyType& key) const {
        return std::hash<KeyType>()(key) % capacity;
    }

public:
    MultiHashTable(size_t cap) : capacity(cap), size(0) {
        table.resize(capacity);
    }

    void insert(const KeyType& key, const ValueType& value) {
        size_t index = hash(key);
        table[index].emplace_back(key, value);
        size++;
    }

    std::vector<ValueType> findAll(const KeyType& key) const {
        std::vector<ValueType> results;
        size_t index = hash(key);
        for (const auto& node : table[index]) {
            if (node.key == key) {
                results.push_back(node.value);
            }
        }
        return results;
    }

    bool contains(const KeyType& key) const {
        size_t index = hash(key);
        for (const auto& node : table[index]) {
            if (node.key == key) {
                return true;
            }
        }
        return false;
    }

    void removeAll(const KeyType& key) {
        size_t index = hash(key);
        auto& bucket = table[index];
        size_t originalSize = bucket.size();

        bucket.remove_if([&key](const HashNode& node) {
            return node.key == key;
            });

        size -= (originalSize - bucket.size());
    }

    void remove(const KeyType& key, const ValueType& value) {
        size_t index = hash(key);
        auto& bucket = table[index];

        for (auto it = bucket.begin(); it != bucket.end(); ) {
            if (it->key == key && it->value == value) {
                it = bucket.erase(it);
                size--;
            }
            else {
                ++it;
            }
        }
    }


    size_t getSize() const {
        return size;
    }

    size_t countValues(const KeyType& key) const {
        size_t count = 0;
        size_t index = hash(key);

        for (const auto& node : table[index]) {
            if (node.key == key) {
                count++;
            }
        }

        return count;
    }

    void clear() {
        table.clear();
        table.resize(capacity);
        size = 0;
    }
};



void FindIS(std::vector<std::vector<uint32_t>>& p20, std::vector<std::vector<uint32_t>>& p17, uint32_t w17, uint32_t w18, uint32_t w23, uint32_t w22) {
    uint32_t rands[8];
    generateRandomUInt32Array(rands, 8);

    uint32_t p18[8];
    uint32_t W17 = w17;
    uint32_t W17_prime;            
    uint32_t W18 = w18;
    uint32_t W18_prime = XOR(w18, w22);
    uint32_t W19;                 
    uint32_t W19_prime;
    uint32_t W21;              
    uint32_t temp;


    for (int i = 0; i < (1 << d_f); ++i) {
        uint32_t ind = (((i & 0b10000) >> 5) << 31) +
            (((i & 0b01100) >> 2) << 22) +
            ((i & 0b00011) << 15);
        W21 = ind;

        
        W17_prime = XOR(W17, W21);
        p18[0] = rands[0] + W17_prime;   
        temp = p18[0] & mask_backward;
        p18[1] = rands[1];              
        p18[1] = p18[1] & 0xe9f5f1f9;
        p18[1] = p18[1] + temp;
        p18[2] = rands[2];              
        p18[3] = rands[3];             
        p18[4] = rands[4] | mask_backward;           
        p18[5] = rands[5];              
        p18[6] = rands[6];              
        p18[7] = rands[7];             

        uint32_t p_tmp[8];
        for (int k = 0; k < 8; k++) p_tmp[k] = p18[k];
        StepFunction(p_tmp, W18, W18_prime, 18);            
        StepFunction(p_tmp, 0xe9f5f1f9, 0xe9f5f1f9, 19);   
        for (int k = 0; k < 8; k++) p20[i][k] = p_tmp[k];
    }

    for (int i = 0; i < (1 << d_b); ++i) {
        uint32_t ind = (((i & 0b1000000000) >> 9) << 28) +
            (((i & 0b0110000000) >> 7) << 25) +
            (((i & 0b0001000000) >> 6) << 19) +
            (((i & 0b0000100000) >> 5) << 17) +
            (((i & 0b0000011100) >> 2) << 9) +
            ((i & 0b0000000011) << 1) + 0xe9f5f1f9;
        W19 = ind;
        W19_prime = XOR(W19, w23);

        p18[0] = W17;                       
        p18[1] = rands[1];                  
        p18[1] = p18[1] & 0xe9f5f1f9;
        p18[1] = p18[1] + temp;                 
        p18[2] = rands[2] - W19_prime;          
        p18[3] = rands[3];                     
        p18[4] = rands[4] | mask_backward;                    
        p18[5] = rands[5];                    
        p18[6] = rands[6] - W19;                
        p18[7] = rands[7];                      

        uint32_t p_tmp[8];
        for (int k = 0; k < 8; k++) p_tmp[k] = p18[k];
        W17_prime = XOR(W17, 0x00000000);
        InvStepFunction(p_tmp, W17, W17_prime, 17);          
        for (int k = 0; k < 8; k++) p17[i][k] = p_tmp[k];
    }
}

void PseudoPreimage_MITM() {



    uint32_t rand15[4];
    generateRandomUInt32Array(rand15, 4);
    uint32_t W17 = rand15[0], W18 = rand15[1], W23 = rand15[2], W22 = rand15[3];

    std::vector<std::vector<uint32_t>> p17((1 << d_b), std::vector<uint32_t>(8, 0));
    std::vector<std::vector<uint32_t>> p20((1 << d_f), std::vector<uint32_t>(8, 0));
    FindIS(p20, p17, W17, W18, W23, W22);


    //Online phase
    uint32_t N_sample = 1 << 14; long long ctr_Pr_b = 0; int ctr_partial_match = 0;
    std::cout << "Number of total samples: " << N_sample << std::endl;

    for (uint32_t sample = 0; sample < N_sample; ++sample)
    {
        uint32_t W[36];
        uint32_t W_prime[36];
        uint32_t rands[6];
        generateRandomUInt32Array(rands, 6);
        W[11] = rands[0], W[14] = rands[1], W[15] = rands[2], W[20] = rands[3];
        W[25] = rands[4];
        W[17] = W17;
        W[18] = W18;
        W[22] = W22;
        W[23] = W23;

        MultiHashTable<uint32_t, int> hashTable((1 << d_f));
        std::vector<std::vector<uint32_t>> auxiTable((1 << d_f), std::vector<uint32_t>(8, 0));
        std::vector<uint32_t> auxiTableW28_prime((1 << d_f), 0);
        std::vector<std::vector<uint32_t>> auxiTableWW((1 << d_f), std::vector<uint32_t>(4, 0));

        //Forward computation
        for (int i = 0; i < (1 << d_f); i++) { //W21
            //message compensation
            uint32_t ind = (((i & 0b10000) >> 5) << 31) +
                (((i & 0b01100) >> 2) << 22) +
                ((i & 0b00011) << 15);
            W[21] = ind;
            W[24] = InvP1(rol(W[21], 15));
            uint32_t W19 = 0xe9f5f1f9;
            uint32_t W16 = InvP1(rol(W19, 7));
            uint32_t W13 = InvP1(rol(W16, 7));
            uint32_t W12 = W19;
            uint32_t W10 = InvP1(rol(W13, 7));
            uint32_t p_tmp[8];
            for (int k = 0; k < 8; k++) p_tmp[k] = p20[i][k];

            // 字扩展
            W[26] = XOR(P1(W10),   XOR(P1(W[17]), XOR(P1(rol(W[23], 15)), XOR(rol(W13, 7),   W[20]))));
            W[27] = XOR(P1(W[11]), XOR(P1(W[18]), XOR(P1(rol(W[24], 15)), XOR(rol(W[14], 7), W[21]))));
            W[28] = XOR(P1(W12),   XOR(P1(W19),   XOR(P1(rol(W[25], 15)), XOR(rol(W[15], 7), W[22]))));
            W[29] = XOR(P1(W13),   XOR(P1(W[20]), XOR(P1(rol(W[26], 15)), XOR(rol(W16, 7),   W[23]))));
            W[30] = XOR(P1(W[14]), XOR(P1(W[21]), XOR(P1(rol(W[27], 15)), XOR(rol(W[17], 7), W[24]))));
            W[31] = XOR(P1(W[15]), XOR(P1(W[22]), XOR(P1(rol(W[28], 15)), XOR(rol(W[18], 7), W[25]))));
            W[32] = XOR(P1(W16),   XOR(P1(W[23]), XOR(P1(rol(W[29], 15)), XOR(rol(W19, 7),   W[26]))));
            W[33] = XOR(P1(W[17]), XOR(P1(W[24]), XOR(P1(rol(W[30], 15)), XOR(rol(W[20], 7), W[27]))));
            W[34] = XOR(P1(W[18]), XOR(P1(W[25]), XOR(P1(rol(W[31], 15)), XOR(rol(W[21], 7), W[28]))));
            W[35] = XOR(P1(W19),   XOR(P1(W[26]), XOR(P1(rol(W[32], 15)), XOR(rol(W[22], 7), W[29]))));



            for (int k = 20; k < 31; k++) {
                W_prime[k] = XOR(W[k], W[k + 4]);
                StepFunction(p_tmp, W[k], W_prime[k], k);
            }

            StepFunctionSP(p_tmp, 0);


            auxiTable[i] = std::vector<uint32_t>(8, 0);
            for (int k = 0; k < 8; k++) auxiTable[i][k] = p_tmp[k]; //p4


            uint32_t A1 = p_tmp[0];
            uint32_t value = A1 & mask_match;
            hashTable.insert(value, i);                 

        }

        for (int i = 0; i < (1 << d_b); ++i) {
            uint32_t ind = (((i & 0b1000000000) >> 9) << 28) +
                (((i & 0b0110000000) >> 7) << 25) +
                (((i & 0b0001000000) >> 6) << 19) +
                (((i & 0b0000100000) >> 5) << 17) +
                (((i & 0b0000011100) >> 2) << 9) +
                ((i & 0b0000000011) << 1) + 0xe9f5f1f9;
            W[19] = ind;
            W[16] = InvP1(rol(W[19], 7));
            W[13] = InvP1(rol(W[16], 7));
            W[12] = W[19];
            W[10] = InvP1(rol(W[13], 7));


            uint32_t p_tmp[8];
            for (int k = 0; k < 8; k++) p_tmp[k] = p17[i][k];
            uint32_t W21 = 0x00000000;
            uint32_t W24 = P1(rol(W21, 15));
        
            W[9] = XOR(InvP1(W[25]), XOR(InvP1(W[19]), XOR(InvP1(rol(W[12], 7)), XOR(W[16], rol(W[22], 15)))));
            W[8] = XOR(InvP1(W[24]), XOR(InvP1(W[18]), XOR(InvP1(rol(W[11], 7)), XOR(W[15], rol(W[21], 15)))));
            W[7] = XOR(InvP1(W[23]), XOR(InvP1(W[17]), XOR(InvP1(rol(W[10], 7)), XOR(W[14], rol(W[20], 15)))));
            W[6] = XOR(InvP1(W[22]), XOR(InvP1(W[16]), XOR(InvP1(rol(W[9], 7)),  XOR(W[13], rol(W[19], 15)))));
            W[5] = XOR(InvP1(W21),   XOR(InvP1(W[15]), XOR(InvP1(rol(W[8], 7)),  XOR(W[12], rol(W[18], 15)))));
            W[4] = XOR(InvP1(W[20]), XOR(InvP1(W[14]), XOR(InvP1(rol(W[7], 7)),  XOR(W[11], rol(W[17], 15)))));
            W[3] = XOR(InvP1(W[19]), XOR(InvP1(W[13]), XOR(InvP1(rol(W[6], 7)),  XOR(W[10], rol(W[16], 15)))));
            W[2] = XOR(InvP1(W[18]), XOR(InvP1(W[12]), XOR(InvP1(rol(W[5], 7)),  XOR(W[9],  rol(W[15], 15)))));
            W[1] = XOR(InvP1(W[17]), XOR(InvP1(W[11]), XOR(InvP1(rol(W[4], 7)),  XOR(W[8],  rol(W[14], 15)))));
            W[0] = XOR(InvP1(W[16]), XOR(InvP1(W[10]), XOR(InvP1(rol(W[3], 7)),  XOR(W[7],  rol(W[13], 15)))));
            for (int j = 16; j >= 0; j--) {
                W_prime[j] = XOR(W[j], W[j + 4]);
            }

            for (int j = 16; j > 5; j--) {
                InvStepFunction(p_tmp, W[j], W_prime[j], j);
            }

            uint32_t pb6[8];
            for (int k = 0; k < 8; k++) pb6[k] = p_tmp[k];

            for (int j = 5; j > 3; j--) {
                InvStepFunction(p_tmp, W[j], W_prime[j], j);
            }


            uint32_t D4 = p_tmp[3];
            uint32_t A1_ = ror(D4, 9);
            uint32_t A1_hat = A1_ - W_prime[0];
            uint32_t value = (A1_hat & mask_match);

            std::vector<int>all_index_W21 = hashTable.findAll(value);
 
            for (int index_W21 : all_index_W21) {
                uint32_t p_tmp_recomputed_0[8];
                for (int k = 0; k < 8; k++) p_tmp_recomputed_0[k] = pb6[k];

                uint32_t W21 = (((index_W21 & 0b10000) >> 5) << 31) +
                    (((index_W21 & 0b01100) >> 2) << 22) +
                    ((index_W21 & 0b00011) << 15);
                //std::cout << W21 << std::endl;
                uint32_t W5 = XOR(InvP1(W21), XOR(InvP1(W[15]), XOR(InvP1(rol(W[8], 7)), XOR(W[12], rol(W[18], 15)))));
                uint32_t W4 = XOR(InvP1(W[20]), XOR(InvP1(W[14]), XOR(InvP1(rol(W[7], 7)), XOR(W[11], rol(W[17], 15)))));
                uint32_t W5_prime = XOR(W5, W[9]);
                uint32_t W4_prime = XOR(W4, W[8]);
                InvStepFunction(p_tmp_recomputed_0, W5, W5_prime, 5); 
                InvStepFunction(p_tmp_recomputed_0, W4, W4_prime, 4); 

                uint32_t D4_recomputed = p_tmp_recomputed_0[3];
                uint32_t A1recomputed_ = ror(D4_recomputed, 9);
                uint32_t A1recomputed_hat = A1recomputed_ - W_prime[0];
                uint32_t value_recomputed = (A1recomputed_hat & mask_match);


                if (value_recomputed == value) {
                    ctr_Pr_b++;

                }
            }
        }
    }
    double pr_f = (double)ctr_Pr_b / (double)(N_sample * (1 << (d_f + d_b - d_m)));
    std::cout << "1. re-estimate Pr_b = Pr[correctly expand three steps in backward]: " << std::fixed << std::setprecision(3) << pr_f << std::endl;
    std::cout << "--------" << ctr_Pr_b << "-----------" << std::endl;

}





int main(int argc, char** argv) {

    PseudoPreimage_MITM();

    return EXIT_SUCCESS;
}
