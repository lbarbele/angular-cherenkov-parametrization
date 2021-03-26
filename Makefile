CXX = g++
CXXFLAGS += `root-config --cflags --libs` -O3 -lgsl

OBJDIR = obj
INCDIR = include
SRCDIR = src

INCLUDES = -I $(INCDIR)
HEADERS  = $(addprefix $(INCDIR)/, cerenkovFitter.h cerenkovHistogramFile.h cerenkovAngularDistribution.h cerenkovCanvas.h cerenkovDeviationProfile.h)

# programs
fit: $(OBJDIR)/fit.o $(OBJDIR)/cerenkovFitter.o $(OBJDIR)/cerenkovHistogramFile.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
	
testParametrization: $(OBJDIR)/testParametrization.o $(OBJDIR)/cerenkovHistogramFile.o $(OBJDIR)/cerenkovAngularDistribution.o $(OBJDIR)/cerenkovFitter.o $(OBJDIR)/cerenkovDeviationProfile.o $(OBJDIR)/cerenkovCanvas.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

# default compilation of cpp files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	@mkdir -p obj
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

# PHONY targets
clean:
	@-rm -fv fit singleFit testParametrization parametrization compareDistributions printParameters
	@-rm -rfv obj

.PHONY: clean
