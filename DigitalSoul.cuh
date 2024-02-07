#include<cstdlib>
#include<stdexcept>
#include<string>
#include<iostream>
#include <cmath>
#include <cstdint>

// Forward declaration of the QN namespace


// Forward declaration of the Complex class template within the QN namespace
namespace QN {
    template <typename T> class Complex;
}

// Forward declaration of the Qudit class template within the QN namespace
namespace QN {
    template <typename T> class Qudit;
}

// Forward declaration of the Gate class template within the QN namespace
namespace QN {
    template <typename T> class Gate;
}

// Forward declaration of the metadata namespace within the QN namespace
namespace QN {
    namespace metadata{};
}

namespace metadata
{
	namespace VHDL
	{
		std::string fixed_point_arithmetic=R"(library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity HalfAdder is Port(a,b:in std_logic; c_out, sum:out std_logic); end HalfAdder;

architecture Behavioral of HalfAdder is
begin
    c_out<=a and b;
    sum<=a xor b;
end Behavioral;
------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity FullAdder is
    Port ( a : in STD_LOGIC;
           b : in STD_LOGIC;
           c_in : in STD_LOGIC;
           c_out : out STD_LOGIC;
           sum : out STD_LOGIC);
end FullAdder;

architecture Behavioral of FullAdder is

begin
    sum<=a xor b xor c_in;
    c_out<=((a xor b) and c_in) or (a and b);
end Behavioral;

------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity FixedPointAdder is -- a+b=c
    generic(N:natural:=32);
    port(a,b:in std_logic_vector(N-1 downto 0); c: out std_logic_vector(N-1 downto 0); overflow:out std_logic);
end FixedPointAdder;

architecture Behavioral of FixedPointAdder is
    signal c_inter: std_logic_vector(N-1 downto 0);
begin
    FA1: entity work.FullAdder port map(a=>a(0), b=>b(0), c_in=>'0', c_out=>c_inter(0), sum=>c(0));
    FA_array: for i in 1 to N-1 generate FA: entity work.FullAdder port map(a=>a(i), b=>b(i),c_in=>c_inter(i-1),c_out=>c_inter(i), sum=>c(i));
    end generate FA_array;
    overflow<=c_inter(N-1);

end architecture Behavioral;
-----------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity FixedPointSubtractor is -- a - b
    generic(N: natural:=32);
    port(a,b:in std_logic_vector(N-1 downto 0); c: out std_logic_vector(N-1 downto 0); overflow: out  std_logic);
end FixedPointSubtractor;

architecture Behavioral of FixedPointSubtractor is
    signal c_inter,b_inter: std_logic_vector(N-1 downto 0);
begin
    b_inter<=not(b);
    FA1: entity work.FullAdder port map(a=>a(0), b=>b_inter(0),c_in=>'1',c_out=>c_inter(0), sum=>c(0));
    FA_array: for i in 1 to N-1 generate FA: entity work.FullAdder port map(a=>a(i), b=>b_inter(i),c_in=>c_inter(i-1),c_out=>c_inter(i), sum=>c(i)); end generate FA_array;
    overflow<=not c_inter(N-1);
end Behavioral;
----------------------------------------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity FixedPointMultiplier_MixedPrecision is
    generic(N:natural:=32);
    port(a,b:in std_logic_vector(N-1 downto 0); c: out std_logic_vector(2*N-1 downto 0));
end FixedPointMultiplier_MixedPrecision;

architecture Behavioral of FixedPointMultiplier_MixedPrecision is
    signal anded: std_logic_vector(N**2-1 downto 0);
    signal sumout: std_logic_vector(N**2+N-1 downto 0);--N rows, N+1 cols
    signal carries: std_logic_vector(N-1 downto 0);
begin
    and_array: for i in 0 to N-1 generate 
        inner_loop: for j in 0 to N-1 generate
            gated: anded(N*i+j)<=a(j) and b(i); 
        end generate inner_loop; 
    end generate and_array;
    carries(0)<='0';
    sumout(N downto 1)<=anded(N-1 downto 0);

    add_array: for i in 1 to N-1 generate
        adder: entity work.FixedPointAdder generic map(N=>N+1) 
        port map(a=>carries(i-1)&sumout((N+1)*(i)-1 downto (N+1)*(i-1)+1),
                                                                        b=>anded(N*(i+1)-1 downto N*i)&'0',
                                                                        overflow=>carries(i),
                                                                        c=>sumout((N+1)*(i+1)-1 downto (N+1)*i) );
    end generate add_array;

    outport1: for i in 0 to N-2 generate
        outport: c(i)<=sumout((N+1)*(i)+1);
    end generate outport1;

    c(2*N-2 downto N-1)<=sumout(N**2+N-1 downto N**2);
    c(2*N-1)<=carries(N-1);
end Behavioral;
----------------------------------------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity FixedPointFMA is
    generic(N:natural:=32);
    port(mul1, mul2, add: in std_logic_vector(N-1 downto 0); fused: out std_logic_vector(N-1 downto 0); overflow: out std_logic);
end FixedPointFMA;


architecture Behavioral of FixedPointFMA is
    signal distillate: std_logic_vector(2*N-1 downto 0);
    signal of_mul, of_add:std_logic;
begin


    multiplier: entity work.FixedPointMultiplier_MixedPrecision generic map(N=>N) port map(a=>mul1, b=>mul2, c=>distillate);
    process(distillate) begin
    if distillate(2*N-1 downto N) = (N-1 downto 0 => '0') then 
        of_mul<='0'; 
    else of_mul <= '1'; 
        end if;
    end process;
    adder: entity work.FixedPointAdder generic map(N=>N) port map(a=>distillate(N-1 downto 0), b=>add, c=>fused, overflow=>of_add);
    
    overflow<=of_mul or of_add;
end architecture Behavioral;)";}
};

namespace kernels
{
    __global__ void LUTx_1(bool *inputs, bool *output, size_t *LogicID, size_t *num_inputs)
    {
     
        size_t recall_index(0);
        for(size_t i(0); i<*num_inputs; i++)
        {
          recall_index |= *(inputs+i)<<i;
        
        };
        *output=((*LogicID)>>recall_index & 1);
    }


    template<typename T, typename s> __global__ void poly(T *a, T operand,T *c,s *N_1)
	{
	c[threadIdx.x]= a[threadIdx.x] * pow(operand, (*N_1)-threadIdx.x);
	}
	
	template<typename T>
	void matmul(T *A, T *B, T *C, size_t N, size_t K, size_t M)
	{
	    for (size_t n = 0; n < N; n++)
	    {
	    	const size_t ind1 = n * K;
		for (size_t m = 0; m < M; m++)
		{    
		    T loader=0;
		    for (size_t k = 0; k < K; ++k)
		    {

		        loader+= A[ind1 + k] * B[k * M + m];
		    }
		    C[n*M+m]=loader;
		}
	    }
	};

}

class LUTx_1 {
private:
    size_t num_inputs;
    size_t imm; 
    size_t LogicID;
    bool residency[4]={true, false, false, false};//0 CPU-C++, 1 NVGPU-CUDA+C++, 2 FPGA-VHDL, 3 QNetwork-C++
    size_t *d_LogicID, *d_num_inputs;

public:
    LUTx_1(size_t inputs, size_t logicID,bool NVGPU=false)
        : num_inputs(inputs), LogicID(logicID)
        {
        	imm=num_inputs-1;
        	
        	if (NVGPU)
        	{
			residency[1]=true;
			cudaMalloc(&d_LogicID, sizeof(size_t));
			cudaMalloc(&d_num_inputs,sizeof(size_t));			
			cudaMemcpy(d_LogicID,&LogicID, sizeof(size_t), cudaMemcpyHostToDevice);
			cudaMemcpy(d_num_inputs,&num_inputs, sizeof(size_t), cudaMemcpyHostToDevice);
        	}
        	

        }

    // Getter and setter methods
    size_t getNumInputs() const { return num_inputs; }
    void setNumInputs(size_t inputs) { num_inputs = inputs; }

    size_t getLogicID() const { return LogicID; }
    void setLogicID(size_t logicID) 
	{ 
		LogicID = logicID;
		if(residency[1])
		{
			cudaMemcpy(d_LogicID,&LogicID, sizeof(size_t), cudaMemcpyHostToDevice);
		}
	}


    //Compute interface

    void computeCPU(bool *inputs, bool *output, uint8_t predecessor_residency=0)
    {
     
        size_t recall_index(0);
        for(size_t i(0); i<num_inputs; i++)
        {
          recall_index |= *(inputs+i)<<i;
        
        };
        *output=(LogicID>>recall_index & 1);
    }
    

   
    void computeGPU(bool *inputs, bool *output, uint8_t predecessor_residency=0)//0: CPU-C++, 1: NVGPU-CUDA+C++, 3:QGraph-C++
    {
	switch(predecessor_residency)
	{
	case 0:
		bool *ibuf;
		cudaMalloc(&ibuf, num_inputs*sizeof(bool));


		
		cudaMemcpy(ibuf, inputs, num_inputs*sizeof(bool), cudaMemcpyHostToDevice);

		
		kernels::LUTx_1<<<1,1>>>(ibuf, output, this->d_LogicID, this->d_num_inputs);
		
		//cudaMemcpy(output, obuf, sizeof(bool),cudaMemcpyDeviceToHost);
		cudaFree(ibuf);
		break;
	case 1:
		kernels::LUTx_1<<<1,1>>>(inputs, output, this->d_LogicID, this->d_num_inputs);
		break;
		
	}
	return;
   
    }
    
    template <typename T=float> QN::Gate<T> UnitaryGen()
    {
    	
    	QN::Complex<T> *ops=new QN::Complex<T>[(2<<this->num_inputs)*(2<<this->num_inputs)];
		size_t maxi=2<<this->num_inputs;
		for(size_t i=0;i<maxi;i++)
		{
            size_t a0=i*maxi;
			for(size_t j=0;j<maxi;j++)
			{
                size_t a1=a0+j;
                ops[a1].re=0;
                ops[a1].im=0;
                
			}
            
		}
	
		size_t maxit=(2<<this->num_inputs)>>1;
		for(size_t i=0;i<maxit;i++)
		{
			//std::cout<<(LogicID>>i & 1)<<'\n';
	        if((LogicID>>i & 1))
	        {
		        ops[((i<<1 )+ 1)*maxi+(i<<1)].re=1.;
		        ops[(i<<1 )*maxi+(i<<1)+1].re=1.;
	        }

	        else
	        {
		        ops[(i<<1)*maxi+(i<<1)].re=1.;
		        ops[((i<<1 )+ 1)*maxi+(i<<1) + 1].re=1.;
	        }	


		}

			return QN::Gate<T>(2<<this->num_inputs,ops,1e8,true);
    };





	std::string entityGen(std::string ename="LUTx_1")
	{
	    std::string a="library IEEE;\nuse IEEE.std_logic_1164.all;\nentity "+ename+" is\nport(\ninputs: in std_logic_vector(0 to "+std::to_string(this->num_inputs-1)+");\no: out\nstd_logic);\nend "+ename+";\narchitecture structural of "+ename+" is\nbegin\nprocess(inputs)\nbegin\ncase inputs is\n";
	    for(size_t i=0;i<(1<<this->num_inputs);i++)
	    {
	        std::string a0="";
	        for(size_t j=0;j<this->num_inputs;j++)
	        {
	            a0+=std::to_string((i>>j)&1);
	        
	        }
	        a+="when \""+a0+"\" => o <= '"+std::to_string(this->LogicID>>i & 1)+"' ;\n";
	    }
	    a+="when others => o <= '0' ;\nend case;\nend process;\nend structural;";
	    
	    return a;
	
	}
    
    
    
    // human interface
    
    std::string ss()
    {
    
      return "LookUp Table|ID="+std::to_string(LogicID)+" with "+std::to_string(num_inputs)+" inputs";
    }
    
    std::string LookUpTable()
	{
        std::string out="LUT";
	out+=std::to_string(this->num_inputs)+" ID="+std::to_string(this->LogicID)+"\n";
	for (size_t i(0);i<this->num_inputs;i++)
	{
		out+="______";
	
	}
	out+="\n";
	for (size_t i(0);i<num_inputs;i++)
	{
		out+="|in"+std::to_string(i);
	
	}
	out+="|out|\n";

	for (size_t i(0);i<(1<<this->num_inputs);i++)
	{
		for (size_t j(this->num_inputs-1);j!=size_t(-1);j--)
		{
			out+="| "+std::to_string(1 & (i>>(imm-j)))+" ";
			
		}
		out+="| "+std::to_string(1 & (this->getLogicID()>>i))+" |\n";
	}
	
	return out;
	
	}



};

template<typename s=size_t,typename T=float> class Polynomial
{
private:
	s N=0, *d_N;
	T* coefficients, *d_coefficients;
	bool residency[4]={true, false, false, false},cached=false;//0 CPU-C++, 1 NVGPU-CUDA+C++, 2 FPGA-VHDL, 3 QNetwork-C++
    T obuf,*d_obuf;
public:
	Polynomial(){};
	//Polynomial(s N, bool NVGPU): N(N) {};
	Polynomial(s N, T* coefficients, bool NVGPU)
	{
		this->N=N; 
		this->coefficients=coefficients;
		if(NVGPU)
		{
			residency[1]=true;

			cudaMalloc(&d_N, sizeof(s));
			cudaMalloc(&d_obuf, sizeof(T));
			cudaMalloc(&(this->d_coefficients), (this->N)*sizeof(T));

			this->N-=1;
			cudaMemcpy(this->d_N,&(this->N),sizeof(s),cudaMemcpyHostToDevice);
			this->N+=1;
			cudaMemcpy(this->d_coefficients, this->coefficients, N*sizeof(T), cudaMemcpyHostToDevice);
			T *hh=new T[this->N];
			cudaMemcpy(hh,d_coefficients,this->N*sizeof(float),cudaMemcpyDeviceToHost );
		}
		
	}

    void setN(s N)
	{
		cached=false;
		this->N=N;
		if(residency[1])
		{
			cudaMemcpy(this->d_N,&(this->N),sizeof(s),cudaMemcpyHostToDevice);
		}
		
	}
    s getN(void) {return N;};
   
    void setCoefficients(T* coefficients)
	{
		cached=false;
		this->coefficients=coefficients;
		if(residency[1])
		{
			cudaMemcpy(this->d_coefficients, this->coefficients, sizeof(T)*(this->N), cudaMemcpyHostToDevice);
		}
	}
    T *getCoefficients(void) {return coefficients;}


    // Compute utils
	void computeCPU(T *ibuf, T *obuf)
	{
		
		T tmp=1;
		for(s i=N-1; i>0; i--)
		{
			(*obuf)+=tmp*coefficients[i];
			tmp*=*ibuf;
		}
		
		*obuf+=tmp*(*coefficients);
	}

	void computeGPU(T *ibuf, T *obuf)
	{
			
		T *ee, *hee;
		cudaMalloc(&ee, N*sizeof(T));

		kernels::poly<T,s><<<1,N>>>(d_coefficients,*ibuf,ee,d_N);
		hee=new T[N];
		cudaDeviceSynchronize();
		cudaMemcpy(hee,ee, N*sizeof(T),cudaMemcpyDeviceToHost);
		for(size_t i=0;i<N;i++)
		{
			(*obuf)+=hee[i];
		}
	}


	// Transpile
	std::string entityGen(std::string name="")
	{
		
		std::string base=R"(library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity CombinationalPolynomialCompute_FixedPoint_)"+name+R"( is
	generic(bit_width:natural:=)"+std::to_string(8*sizeof(T))+R"(; order:natural:=)"+std::to_string(N-1)+R"();
	Port(x:in std_logic_vector(bit_width-1 downto 0); y:out std_logic_vector(bit_width-1 downto 0); overflow: out std_logic);
end CombinationalPolynomialCompute_FixedPoint_)"+name+R"(;

architecture Behavioral of CombinationalPolynomialCompute_FixedPoint_)"+name+R"( is
	signal pre_multiplied, pre_fma, coeff: std_logic_vector(bit_width*(order+1)-1 downto 0);
	signal upstreamed:  std_logic_vector(2*bit_width*(order)-1 downto 0);
	signal fma_of: std_logic_vector(order-1 downto 0);

begin)";
		std::string aa;
		
		for(size_t i=0;i<N;i++)
		{
			aa=std::string(8*sizeof(T), '0');
			for(size_t j=0;j<8*sizeof(T);j++)
			{
				if(coefficients[N-1-i]&(1<<j)){aa[8*sizeof(T)-j-1]='1';}
				
			}
			base+="\n\tcoeff("+std::to_string((i+1)*8*sizeof(T)-1)+" downto "+std::to_string(i*8*sizeof(T))+")<=\""+aa+'"'+';';
			
		}
		base+=R"(
pre_fma(bit_width-1 downto 0)<=(bit_width-1 downto 0 => '0');
pre_multiplied(bit_width-1 downto 0)<=(bit_width-2 downto 0 => '0')&'1';
	multiply_next: for i in 1 to order generate
		aa:entity work.FixedPointMultiplier_MixedPrecision 
			generic map(N=>bit_width)
			port map(a=>x,
					b=>pre_multiplied((i)*bit_width-1 downto (i-1)*bit_width),
					c=>upstreamed(2*bit_width*(i)-1 downto 2*bit_width*(i-1)));
		pre_multiplied((i+1)*bit_width-1 downto (i)*bit_width)<=upstreamed(2*bit_width*(i)-1-bit_width downto 2*bit_width*(i-1));
	end generate multiply_next; 

fma_array: for i in 1 to order generate
	fma: entity work.FixedPointFMA generic map(N=>bit_width) 
	   port map(overflow=>fma_of(i-1),
	            add=>pre_fma(i*bit_width-1 downto (i-1)*bit_width),
	            mul1=>pre_multiplied((i)*bit_width-1 downto (i-1)*bit_width),
	            mul2=>coeff(i*bit_width-1 downto (i-1)*bit_width),
	            fused=>pre_fma((i+1)*bit_width-1 downto (i)*bit_width));
end generate fma_array;

process(fma_of) begin
    if fma_of/=(order-1 downto 0 => '0') then overflow<='1'; else overflow<='0'; end if;
end process;

y<=pre_fma(bit_width*(order+1)-1 downto bit_width*order);
	
end Behavioral;
	)";
	return base;
	}

	// Human Interface
	std::string ss()
	{
		std::string a=std::to_string(N-1);
		switch((N-1)%10)
		{
			case 1:
				a+="st"; break;
			case 2:
				a+="nd"; break;
			case 3:
				a+="rd"; break;
			default:
				a+="th";
		}

		a+=" order Polynomial P(x)=";
		for(s i=N-1;i>0;i--)
		{
			a+=std::to_string(*(coefficients+this->N-i-1))+'*'+'x'+'^'+std::to_string(i)+'+';
		}
        a+=std::to_string(coefficients[N-1]);

		return a;
	};

    Polynomial operator+ (Polynomial<s,T>& operand)
    {
        bool i_am_max_n=(this->N >= operand.N)?true:false;
        s maxn;
        if (i_am_max_n)
        {
            maxn=this->N;
            T *new_coeff=new T[this->N];
            for(s i=0;i<operand.N;i++)
            {
                *(new_coeff+i)=*(this->coefficients + i) + *(operand.coefficients + i);
            }
            for(s i=operand.N;i<this->N;i++)
            {
                *(new_coeff+i)=*(this->coefficients+i);
            }
            return Polynomial<s,T>(maxn, new_coeff, false);
        }
        else
        {
            T *new_coeff=new T[operand.N];
            maxn=operand.N;
            for(s i=0;i<this->N;i++)
            {
                *(new_coeff+i)=*(this->coefficients + i) + *(operand.coefficients + i);
            }
            for(s i=this->N;i<operand.N;i++)
            {
                *(new_coeff+i)=*(operand.coefficients+i);
            }
            return Polynomial<s,T>(maxn, new_coeff, false);
        }   
    };
    
    Polynomial operator- (Polynomial<s,T>& operand)
    {
        bool i_am_max_n=(this->N >= operand.N)?true:false;
        s maxn;
        if (i_am_max_n)
        {
            maxn=this->N;
            T *new_coeff=new T[this->N];
            for(s i=0;i<operand.N;i++)
            {
                *(new_coeff+i)=*(this->coefficients + i) - *(operand.coefficients + i);
            }
            for(s i=operand.N;i<this->N;i++)
            {
                *(new_coeff+i)=*(this->coefficients+i);
            }
            return Polynomial<s,T>(maxn, new_coeff, false);
        }
        else
        {
            T *new_coeff=new T[operand.N];
            maxn=operand.N;
            for(s i=0;i<this->N;i++)
            {
                *(new_coeff+i)=*(this->coefficients + i) - *(operand.coefficients + i);
            }
            for(s i=this->N;i<operand.N;i++)
            {
                *(new_coeff+i)=-*(operand.coefficients+i);
            }
            return Polynomial<s,T>(maxn, new_coeff, false);
        }   
    };

    Polynomial<s,T> operator*(Polynomial<s,T>& operand)
    {
        s maxn=(this->N)+(operand.N)-1;
        T *new_coeff=new T[maxn];
        for(s i=0; i<this->N;i++)
        {
            for(s j=0; j<operand.N;j++)
            {

                new_coeff[i+j]+=(this->coefficients[i])*(operand.coefficients[j]);
                
            }
            
        }

        return Polynomial<s,T>(maxn, new_coeff,false);
    }

    Polynomial<s,T> operator/(Polynomial<s,T>& operand)
    {
        s maxn=(this->N)-(operand.N)+1;
        T *new_coeff=new T[maxn];
        T *tmp_coeff=new T[this->N];

        for(size_t k=0;k<this->N;k++)
        {
            tmp_coeff[k]=this->coefficients[k];

        }

        for(size_t i=0;i<maxn;i++)
        {
            
            T mul=(tmp_coeff[i])/(*(operand.coefficients));
            new_coeff[i]=mul;
            for(size_t j=0;j<operand.N;j++)
            {
                tmp_coeff[i+j]-=mul*operand.coefficients[j];
            }

            
        }    
        return Polynomial(maxn,new_coeff, false);

        }

};


namespace QN
{
	template<typename T=float> class Complex
	{
	public:
		T re, im;
		Complex(void){}; // Default constructor
		Complex(T re, T im){this->re=re; this->im=im;};
		
		void operator() (T re, T im){this->re=re; this->im=im;};
		
		Complex<T> operator +(const Complex &B)
		{
			return Complex<T>(this->re+B.re, this->im+B.im);
		};
		void operator +=(const Complex &B)
		{
			this->re+=B.re; this->im+=B.im;
		};

		Complex<T> operator -(const Complex &B)
		{
			return Complex<T>(this->re-B.re, this->im-B.im);
		};
		void operator -=(const Complex &B)
		{
			this->re-=B.re; this->im-=B.im;
		}
		
		Complex<T> operator*(const Complex &B)
		{
			return Complex<T>(this->re*B.re-this->im*B.im, this->re*B.im+this->im*B.re);
		};
		void operator*=(const Complex &B)
		{
			 this->re=this->re*B.re-this->im*B.im;
			 this->im=this->re*B.im+this->im*B.re;
		};
		
		Complex<T> operator/(const Complex &B)
		{
			T magB=B.re*B.re + B.im*B.im;
			if(magB!=0){return Complex<T>((this->re*B.re + this->im * B.im)/magB, (this->im * B.re - this->re * B.im)/magB);}
			else{throw std::runtime_error("ZeroDivisionError: Invaild value encountered in Complex division");}
		};
		void operator/=(const Complex &B)
		{
			T magB=B.re*B.re + B.im*B.im;
			if(magB!=0){this->re=(this->re*B.re + this->im * B.im)/magB; this->im=(this->im * B.re - this->re * B.im)/magB;}
			else{throw std::runtime_error("ZeroDivisionError: Invaild value encountered in Complex division");}
		};
		
		T magnitude(void)
		{
			
			return sqrt(this->re*this->re+this->im*this->im);
		};
		
		T arg(void)
		{
			return atan2(this->im,this->re);
		};
		
		Complex<T> conj()
		{
			return Complex<T>(this->re, -1*(this->im));
		}
		
		std::string ss(void)
		{
			if(this->im>=0)
			{
				return std::to_string(this->re)+"+"+std::to_string(this->im)+"j";
			}
			else
			{
				return std::to_string(this->re)+"-"+std::to_string(-1*this->im)+"j";
			};
		};
	};



	template <typename T=float> class Qudit
	{

		size_t N_system=2;
		Complex<T> *statevector=nullptr;
	public:
		Qudit(){};
		Qudit(size_t N_system, bool allocate_now=true)
		{
			this->N_system=N_system;
			if(allocate_now){this->statevector=new Complex<T>[this->N_system];this->statevector[0](1,0);};
			
		};
		
		Qudit(size_t N_system, Complex<T> *new_statevector, T epsilon=1e-8, bool trust=false)
		{
			if(trust){delete [] this->statevector;this->statevector=new_statevector;}
			else
				{
				T cum_sum=0;
				for(size_t i=0;i<this->N_system;i++)
				{
					cum_sum+=new_statevector[i].magnitude();
				}
			
				if(abs(cum_sum-1)>epsilon)
				{
					throw std::runtime_error("DenormalizedPureState: Invalid value encountered in new_statevector. The new_statevector must be normalized to one. Consider to normalize or modify components");
				}
			
				delete [] this->statevector;
				this->statevector=new_statevector;
				}
		};
		
		void oneHot(size_t N)
		{
			for(size_t i=0;i<this->N_system;i++)
			{
				this->statevector[i](0,0);
			}
			this->statevector[N](1,0);
		};
		
		void loadStatevector(Complex<T> *new_statevector, T epsilon=1e-8, bool trust=false)
		{
			if(trust){delete [] this->statevector;this->statevector=new_statevector;}
			else
				{
				T cum_sum=0;
				for(size_t i=0;i<this->N_system;i++)
				{
					cum_sum+=new_statevector[i].magnitude();
				}
			
				if(abs(cum_sum-1)>epsilon)
				{
					throw std::runtime_error("DenormalizedPureState: Invalid value encountered in new_statevector. The new_statevector must be normalized to one. Consider to normalize or modify components");
				}
			
				delete [] this->statevector;
				this->statevector=new_statevector;
				}
		};
		
		void freeStatevector(void)
		{	
			delete [] (this->statevector);
			this->statevector=nullptr;
		};
		
		Complex<T> *Psi(void)
		{return (this->statevector);};
		
		size_t numStates(void)
		{
			return this->N_system;
		};
		
		std::string ss(bool verbosity=false)
		{
			if(verbosity)
			{
				std::string rs=std::to_string(this->N_system)+"-level Qudit\n\u03A8 = (";
				for(size_t i=0;i<this->N_system-1;i++)
				{
					rs+=this->statevector[i].ss()+", ";
				}
				rs+=this->statevector[this->N_system-1].ss()+") ";
				return rs;
			}
			
			else{return std::to_string(this->N_system)+"-level Qudit";};
		};
	};

	template<typename T> bool isUnitary(Complex<T> *A, size_t &operandDim, T epsilon=1e-8)
	{
		
		for(size_t i=0;i<operandDim;i++)
		{
			for(size_t j=0;j<operandDim; j++)
			{
				Complex<T> current(0,0);
				for(size_t k=0;k<operandDim;k++)
				{
					current+=A[i*operandDim+k]*(A[j*operandDim+k].conj());
				}//j,k -> kj
				if(i==j){if(abs(current.re-1)<epsilon && abs(current.im)<epsilon){} else{return false;}}
				else{if(abs(current.re)<epsilon && abs(current.im)<epsilon){}else{return false;}}		
			}
		}
		
		return true;
	}

	template<typename T=float> void sqMatmul(Complex<T> *U1, Complex<T> *U2, Complex<T> *U3, size_t operandDim)
	{
		for(size_t i=0; i<operandDim;i++)
		{
			size_t ai=i*operandDim;
			for(size_t j=0;j<operandDim;j++)
			{
				for(size_t k=0;k<operandDim;k++)
				{
					U3[ai+j]+=U1[ai+k]*U2[k*operandDim+j];
				}
				
			}
		}
		
	}
	template<typename T=float> void sqMVDot(Complex<T> *U, Complex<T> *S1, Complex<T> *S2, size_t operandDim)
	{
		for(size_t i=0;i<operandDim;i++)
		{
			size_t ai=operandDim*i;
			for(size_t j=0;j<operandDim;j++)
			{
				S2[i]+=U[ai+j]*S1[j];
			}
			
		}
	}

	template <typename T=float> class Gate
	{
	private:
		Complex<T> *data=nullptr;
		size_t operandDim=2;
	public:
		Gate(){};
		Gate(size_t operandDim){this->operandDim=operandDim;};
		Gate(size_t OperandDim, Complex<T> *Operator, T epsilon=1e-8, bool trust=false)
		{
			this->operandDim=OperandDim;
			if(trust)
			{
				this->data=Operator;
			}
			else
			{
				if(isUnitary(Operator,this->operandDim, epsilon)){this->data=Operator;}		
				else{throw std::runtime_error("Operator is not Unitary. Try to change either coefficients or accomodate a change for epsilon so that it may not be detected precisely");}
				
			}
		};

		void loadOperator(Complex<T> *Operator, T epsilon=1e-8, bool trust=false)
		{
			if(trust)
			{
				this->data=Operator;
			}
			else
			{
				if(isUnitary(Operator,this->operandDim, epsilon)){this->data=Operator;}		
				else{throw std::runtime_error("Operator is not Unitary. Try to change either coefficients or epsilon so that it may not be detected precisely");}
				
			}
		};
		void transform(Qudit<T> &sv, bool trust=false)// CONTINUE	FLAG
		{
			if(!trust)
			{
				if(sv.numStates()!=this->operandDim)
				{
					throw std::runtime_error("Unmatched dimensions of operator and operand. Gate is "+std::to_string(this->operandDim)+"-level and Qudit is "+std::to_string(sv.numStates())+"-level system");
				}
			}
			
			Complex<T> *s1=new Complex<T>[this->operandDim];
			Complex<T> *s2=sv.Psi();
			for(size_t i=0;i<this->operandDim;i++)
			{
				s1[i]=s2[i];
			}
			
			for(size_t i=0;i<this->operandDim;i++)
			{
				s2[i].re=0.;
				s2[i].im=0.;
			}
		
			sqMVDot(this->data, s1, s2, this->operandDim);
		}
		
		std::string ss()
		{
		  std::string out="Quantum Gate on "+std::to_string(this->operandDim) +"-level system \nU=[";
		  for(size_t i=0;i<this->operandDim;i++)
		  {
		    out+='[';
		    for(size_t j=0;j<this->operandDim-1;j++)
		      {
			out+=this->data[i*this->operandDim+j].ss()+',';
		      
		      }
		      out+=this->data[(i+1)*this->operandDim-1].ss()+']'+','+'\n';
		  
		  }
		  out+=']';
		return out;
		};
		
	};
	
	namespace metadata
	{
		typedef float native;
		Complex<native> Xarray[4]={Complex<native>(0,0),Complex<native>(1,0), Complex<native>(1,0), Complex<native>(0,0) } ;

		Complex<native> Yarray[4]={Complex<native>(0,0),Complex<native>(0,-1), Complex<native>(0,1), Complex<native>(0,0) } ;

		Complex<native> Zarray[4]={Complex<native>(1,0),Complex<native>(0,0), Complex<native>(0,0), Complex<native>(-1,0) } ;
	}

};
