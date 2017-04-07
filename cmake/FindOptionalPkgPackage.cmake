# - FindOptionalPkgPackage
# Enable or disable optional packages. Also force optional packages.
#
#  This module defines the following variables:
#
#
#
#  
# Copyright (c) 2017  Frank Uhlig (uhlig.frank@gmail.com)
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 3 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.  
#  
# This file incorporates work covered by the following copyright and  
# permission notice:  
#
#=============================================================================
# Copyright 2011 Nils Andresen <nils@nils-andresen.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================


macro(find_optional_pkg_package _normal_package)
	STRING(TOUPPER ${_normal_package} _upper_package)
    STRING(TOLOWER ${_normal_package} _lower_package)

	OPTION(ENABLE_${_upper_package} "Force dependencies to ${_normal_package}" OFF)
	OPTION(DISABLE_${_upper_package} "Never depend on ${_normal_package}" OFF)

	if(ENABLE_${_upper_package})
		pkg_search_module(${_normal_package} ${_lower_package} REQUIRED)
	elseif(NOT DISABLE_${_upper_package})
		pkg_search_module(${_normal_package} ${_lower_package})
	endif(ENABLE_${_upper_package})
endmacro(find_optional_pkg_package)
