CXX ?= g++

CXXFLAGS := -std=c++17 -O3 -w
INCLUDES := -Iinfrastructure/include -Iksecforest/include

BUILD_DIR := build

COMMON_SRC := \
	ksecforest/src/ksecforest_mba.cpp \
	infrastructure/src/KDeltaTrussDecomp.cpp

TEST_SRCS := \
	ksecforest/test/constructionTest.cpp \
	ksecforest/test/addEdgeTest.cpp \
	ksecforest/test/queryTest.cpp

TESTS := constructionTest addEdgeTest queryTest

COMMON_OBJ := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(COMMON_SRC))
TEST_OBJ   := $(patsubst ksecforest/test/%.cpp,$(BUILD_DIR)/test/%.o,$(TEST_SRCS))

DEPS := $(COMMON_OBJ:.o=.d) $(TEST_OBJ:.o=.d)

# ======================
# color output
# ======================
GREEN=\033[0;32m
BLUE=\033[0;34m
YELLOW=\033[0;33m
RED=\033[0;31m
NC=\033[0m

# ======================
# silent build (no g++ logs)
# ======================
.SILENT:

.PHONY: all clean

# ======================
# default
# ======================
all: banner $(TESTS)
	@echo "$(GREEN)[DONE] All targets built successfully ✔$(NC)"

# ======================
# banner
# ======================
banner:
	@echo "$(BLUE)======================================$(NC)"
	@echo "$(BLUE)      KSECFOREST BUILD SYSTEM        $(NC)"
	@echo "$(BLUE)======================================$(NC)"

# ======================
# link rules (FINAL OUTPUT ONLY)
# ======================
constructionTest: $(BUILD_DIR)/test/constructionTest.o $(COMMON_OBJ)
	@echo "$(YELLOW)[LINK] constructionTest$(NC)"
	@$(CXX) $^ -o ksecforest/test/$@.exe
	@echo "$(GREEN)[OK] constructionTest ✔$(NC)"

addEdgeTest: $(BUILD_DIR)/test/addEdgeTest.o $(COMMON_OBJ)
	@echo "$(YELLOW)[LINK] addEdgeTest$(NC)"
	@$(CXX) $^ -o ksecforest/test/$@.exe
	@echo "$(GREEN)[OK] addEdgeTest ✔$(NC)"

queryTest: $(BUILD_DIR)/test/queryTest.o $(COMMON_OBJ)
	@echo "$(YELLOW)[LINK] queryTest$(NC)"
	@$(CXX) $^ -o ksecforest/test/$@.exe
	@echo "$(GREEN)[OK] queryTest ✔$(NC)"

# ======================
# compile COMMON
# ======================
$(BUILD_DIR)/ksecforest/src/%.o: ksecforest/src/%.cpp
	@mkdir -p $(dir $@)
	@echo "$(BLUE)[CC] $<$(NC)"
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

$(BUILD_DIR)/infrastructure/src/%.o: infrastructure/src/%.cpp
	@mkdir -p $(dir $@)
	@echo "$(BLUE)[CC] $<$(NC)"
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

# ======================
# compile TEST
# ======================
$(BUILD_DIR)/test/%.o: ksecforest/test/%.cpp
	@mkdir -p $(dir $@)
	@echo "$(BLUE)[CC] $<$(NC)"
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

# ======================
# auto deps
# ======================
-include $(DEPS)

# ======================
# clean
# ======================
clean:
	@echo "$(RED)[CLEAN] removing build artifacts$(NC)"
	@rm -rf $(BUILD_DIR) ksecforest/test/*.exe