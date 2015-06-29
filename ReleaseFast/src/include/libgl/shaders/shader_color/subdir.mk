################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/include/libgl/shaders/shader_color/CShaderColor.cpp 

OBJS += \
./src/include/libgl/shaders/shader_color/CShaderColor.o 

CPP_DEPS += \
./src/include/libgl/shaders/shader_color/CShaderColor.d 


# Each subdirectory must supply rules for building sources it contributes
src/include/libgl/shaders/shader_color/%.o: ../src/include/libgl/shaders/shader_color/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -I"/home/martin/workspace/sweet/src/include" -O3 -g3 -pg -Wall `pkg-config --cflags sdl2` `pkg-config freetype2 --cflags` -fopenmp -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


