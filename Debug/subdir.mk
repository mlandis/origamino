################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Alignment.cpp \
../Code.cpp \
../CondLikes.cpp \
../MbBitfield.cpp \
../MbEigensystem.cpp \
../MbMath.cpp \
../MbRandom.cpp \
../MbTransitionMatrix.cpp \
../Mcmc.cpp \
../Menu.cpp \
../ModelPath.cpp \
../ModelSeq.cpp \
../Parm.cpp \
../Parm_context.cpp \
../Parm_freqs.cpp \
../Parm_kappa.cpp \
../Parm_omega.cpp \
../Parm_pairwise.cpp \
../Parm_scaling.cpp \
../Parm_solv.cpp \
../Parm_tree.cpp \
../Path.cpp \
../Settings.cpp \
../TiProbs.cpp \
../main.cpp 

OBJS += \
./Alignment.o \
./Code.o \
./CondLikes.o \
./MbBitfield.o \
./MbEigensystem.o \
./MbMath.o \
./MbRandom.o \
./MbTransitionMatrix.o \
./Mcmc.o \
./Menu.o \
./ModelPath.o \
./ModelSeq.o \
./Parm.o \
./Parm_context.o \
./Parm_freqs.o \
./Parm_kappa.o \
./Parm_omega.o \
./Parm_pairwise.o \
./Parm_scaling.o \
./Parm_solv.o \
./Parm_tree.o \
./Path.o \
./Settings.o \
./TiProbs.o \
./main.o 

CPP_DEPS += \
./Alignment.d \
./Code.d \
./CondLikes.d \
./MbBitfield.d \
./MbEigensystem.d \
./MbMath.d \
./MbRandom.d \
./MbTransitionMatrix.d \
./Mcmc.d \
./Menu.d \
./ModelPath.d \
./ModelSeq.d \
./Parm.d \
./Parm_context.d \
./Parm_freqs.d \
./Parm_kappa.d \
./Parm_omega.d \
./Parm_pairwise.d \
./Parm_scaling.d \
./Parm_solv.d \
./Parm_tree.d \
./Path.d \
./Settings.d \
./TiProbs.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


