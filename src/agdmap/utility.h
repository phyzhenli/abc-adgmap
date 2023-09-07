/**CFile****************************************************************

  FileName    [utility.h]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Agdmap.]

  Synopsis    [External declarations.]

  Author      [Longfei Fan]

  Affiliation [Fudan University]

  Date        [August 1, 2023]

  Vision      [ 1.0 ]

***********************************************************************/
#pragma once

#include <vector>
#include <list>
#include <queue>
#include <map>
#include <deque>
#include <stack>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>
#include <functional>
#include <future>
#include <utility>
#include <cassert>
#include <ctime>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <memory>
#include <algorithm>
#include <atomic>
#include <unistd.h>

namespace utility
{

  template <class TYPE>
  class QMP
  {
  public:
    struct Implicant
    {
      int implicant;
      std::string minterms;
      std::vector<int> mints;
      int mask;
      std::string bits;
      int ones;
      bool used;
      Implicant(int var_cnt,
                int i = 0, std::vector<int> min = std::vector<int>(),
                std::string t = "",
                int m = 0, bool u = false) : implicant(i), mask(m), ones(0), used(u)
      {
        if (t == "")
        {
          std::stringstream ss;
          ss << 'm' << i;
          minterms = ss.str();
        }
        else
          minterms = t;
        if (min.empty())
          mints.push_back(i);
        else
          mints = min;
        int bit = 1 << var_cnt;
        while (bit >>= 1)
          if (m & bit)
            bits += '-';
          else if (i & bit)
          {
            bits += '1';
            ++ones;
          }
          else
            bits += '0';
      }
      bool operator<(const Implicant &b) const
      {
        return ones < b.ones;
      }
      std::vector<int> cat(const Implicant &b)
      {
        std::vector<int> v = mints;
        v.insert(v.end(), b.mints.begin(), b.mints.end());
        return v;
      }
      class QMP;
    };

  public:
    QMP() : output_("F") {}

    void clear() { vars_.clear(); }
    void addVar(const std::string &var_name) { vars_.push_back(var_name); }
    void set_output(const std::string &name) { output_ = name; }

    int var_cnt() { return static_cast<int>(vars_.size()); }
    std::string getMinExpression(TYPE true_value)
    {
      const float b = 2.0f;
      const float a = sizeof(TYPE);
      int max_var_cnt = static_cast<int>(floor(log(a) / log(b))) + 3;
      // int max_var_cnt = static_cast<int>(
      //     floor(log(sizeof(unsigned int)) / log(2))) + 3;
      if (var_cnt() > max_var_cnt)
        return "";

      int combs = 1 << var_cnt();
      std::vector<int> minterms;
      std::vector<Implicant> implicants;
      for (int mint = 0; mint < combs; mint++)
      {
        if ((true_value >> mint) & 1)
        {
          implicants.push_back(Implicant(var_cnt(), mint));
          minterms.push_back(mint);
        }
      }
      if (!minterms.size())
      {
        std::string str_exp = output_ + "=0";
        return str_exp;
      }
      std::sort(minterms.begin(), minterms.end());
      minterms.erase(std::unique(minterms.begin(), minterms.end()), minterms.end());
      std::sort(implicants.begin(), implicants.end());
      // print(implicants);

      std::vector<Implicant> aux;
      std::vector<Implicant> primes;
      while (implicants.size() > 1)
      {
        for (size_t i = 0; i < implicants.size() - 1; ++i)
          for (size_t j = implicants.size() - 1; j > i; --j)
            if (implicants[j].bits == implicants[i].bits)
              implicants.erase(implicants.begin() + j);
        aux.clear();
        for (size_t i = 0; i < implicants.size() - 1; ++i)
          for (size_t j = i + 1; j < implicants.size(); ++j)
            if (implicants[j].ones == implicants[i].ones + 1 && implicants[j].mask == implicants[i].mask && count1s(implicants[i].implicant ^ implicants[j].implicant) == 1)
            {
              implicants[i].used = true;
              implicants[j].used = true;
              aux.push_back(Implicant(var_cnt(),
                                      implicants[i].implicant,
                                      implicants[i].cat(implicants[j]),
                                      implicants[i].minterms + ',' + implicants[j].minterms,
                                      (implicants[i].implicant ^ implicants[j].implicant) | implicants[i].mask));
            }
        for (size_t i = 0; i < implicants.size(); ++i)
          if (!implicants[i].used)
            primes.push_back(implicants[i]);
        implicants = aux;
        std::sort(implicants.begin(), implicants.end());
        // print(implicants);
      }
      for (size_t i = 0; i < implicants.size(); ++i)
        primes.push_back(implicants[i]);
      if (primes.back().mask == combs - 1)
      {
        std::string str_exp = output_ + "=1";
        return str_exp;
      }

      std::vector<std::vector<bool>> table(primes.size(),
                                           std::vector<bool>(minterms.size()));
      // bool table[primes.size()][minterms.size()];
      for (size_t i = 0; i < primes.size(); ++i)
        for (size_t k = 0; k < minterms.size(); ++k)
          table[i][k] = false;
      for (size_t i = 0; i < primes.size(); ++i)
        for (size_t j = 0; j < primes[i].mints.size(); ++j)
          for (size_t k = 0; k < minterms.size(); ++k)
            if (primes[i].mints[j] == minterms[k])
              table[i][k] = true;

      // Petrick
      assert(primes.size() < 64);
      std::vector<uint64_t> M0, M1;
      for (size_t i = 0; i < primes.size(); ++i)
        if (table[i][0])
          M0.push_back(1ull << i);
      for (size_t k = 1; k < minterms.size(); ++k)
      {
        M1.clear();
        for (size_t i = 0; i < primes.size(); ++i)
          if (table[i][k])
            M1.push_back(1ull << i);
        mul(M0, M1);
      }
      int min = count1s(M0[0]);
      size_t ind = 0;
      for (size_t i = 1; i < M0.size(); ++i)
      {
        if (min > count1s(M0[i]))
        {
          min = count1s(M0[i]);
          ind = i;
        }
      }

      /*bool f = false;
      std::string str_exp;
      for (size_t j = 0; j < M0.size(); ++j) {
        str_exp = output_ + "=";
        f = false;
        for (size_t i = 0; i < primes.size(); ++i) {
          if (M0[j] & (1 << i)) {
            if (f)
              str_exp += "+";
            f = true;
            str_exp += getName(primes[i]);
          }
        }
      }*/

      // minimal solution
      bool f = false;
      std::string str_exp = output_ + "=";
      for (size_t i = 0; i < primes.size(); ++i)
      {
        if (M0[ind] & (1ull << i))
        {
          if (f)
            str_exp += "+";
          f = true;
          // str_exp += getName(primes[i]);
          str_exp += addBracket(getName(primes[i]));
        }
      }
      return str_exp;
    }

    std::string addBracket(std::string exp)
    {
      std::string out_str = "";
      for (size_t i = 0; i < exp.length();)
      {
        std::string temp_name = "";
        char token_1 = '\0';
        char token_2 = '\0';
        if (exp[i] == '*')
          token_1 = exp[i++];
        if (exp[i] == '~')
          token_2 = exp[i++];
        bool find = false;
        while (i < exp.length() && find == false)
        {
          temp_name += exp[i++];
          for (int j = 0; j < var_cnt(); ++j)
          {
            if (temp_name == vars_[j])
            {
              find = true;
              break;
            }
          }
        }
        if (find == true)
        {
          if (token_2 == '~')
            temp_name = "(~" + temp_name + ")";
          if (token_1 == '*')
            out_str = "(" + out_str + "*" + temp_name + ")";
          else
            out_str = temp_name;
        }
        else
        {
          assert(0);
          return "";
        }
      }
      return out_str;
    }

    std::vector<bool> getUsedVars(TYPE true_value)
    {
      const float b = 2.0f;
      const float a = sizeof(TYPE);
      int max_var_cnt = static_cast<int>(floor(log(a) / log(b))) + 3;
      // int max_var_cnt = static_cast<int>(
      //     floor(log(sizeof(unsigned int)) / log(2))) + 3;
      std::vector<bool> used_var_flags(var_cnt(), false);
      if (var_cnt() > max_var_cnt)
        return used_var_flags;

      int combs = 1 << var_cnt();
      std::vector<int> minterms;
      std::vector<Implicant> implicants;
      for (int mint = 0; mint < combs; mint++)
      {
        if ((true_value >> mint) & 1)
        {
          implicants.push_back(Implicant(var_cnt(), mint));
          minterms.push_back(mint);
        }
      }
      if (!minterms.size())
      {
        return used_var_flags;
      }
      std::sort(minterms.begin(), minterms.end());
      minterms.erase(std::unique(minterms.begin(), minterms.end()), minterms.end());
      std::sort(implicants.begin(), implicants.end());
      // print(implicants);

      std::vector<Implicant> aux;
      std::vector<Implicant> primes;
      while (implicants.size() > 1)
      {
        for (size_t i = 0; i < implicants.size() - 1; ++i)
          for (size_t j = implicants.size() - 1; j > i; --j)
            if (implicants[j].bits == implicants[i].bits)
              implicants.erase(implicants.begin() + j);
        aux.clear();
        for (size_t i = 0; i < implicants.size() - 1; ++i)
          for (size_t j = i + 1; j < implicants.size(); ++j)
            if (implicants[j].ones == implicants[i].ones + 1 && implicants[j].mask == implicants[i].mask && count1s(implicants[i].implicant ^ implicants[j].implicant) == 1)
            {
              implicants[i].used = true;
              implicants[j].used = true;
              aux.push_back(Implicant(var_cnt(),
                                      implicants[i].implicant,
                                      implicants[i].cat(implicants[j]),
                                      implicants[i].minterms + ',' + implicants[j].minterms,
                                      (implicants[i].implicant ^ implicants[j].implicant) | implicants[i].mask));
            }
        for (size_t i = 0; i < implicants.size(); ++i)
          if (!implicants[i].used)
            primes.push_back(implicants[i]);
        implicants = aux;
        std::sort(implicants.begin(), implicants.end());
        // print(implicants);
      }
      for (size_t i = 0; i < implicants.size(); ++i)
        primes.push_back(implicants[i]);
      if (primes.back().mask == combs - 1)
      {
        return used_var_flags;
      }

      std::vector<std::vector<bool>> table(primes.size(),
                                           std::vector<bool>(minterms.size()));
      // bool table[primes.size()][minterms.size()];
      for (size_t i = 0; i < primes.size(); ++i)
        for (size_t k = 0; k < minterms.size(); ++k)
          table[i][k] = false;
      for (size_t i = 0; i < primes.size(); ++i)
        for (size_t j = 0; j < primes[i].mints.size(); ++j)
          for (size_t k = 0; k < minterms.size(); ++k)
            if (primes[i].mints[j] == minterms[k])
              table[i][k] = true;

      // Petrick
      std::vector<size_t> M0, M1;
      for (size_t i = 0; i < primes.size(); ++i)
        if (table[i][0])
          M0.push_back(1 << i);
      for (size_t k = 1; k < minterms.size(); ++k)
      {
        M1.clear();
        for (size_t i = 0; i < primes.size(); ++i)
          if (table[i][k])
            M1.push_back(1 << i);
        mul(M0, M1);
      }
      int min = count1s(M0[0]);
      size_t ind = 0;
      for (size_t i = 1; i < M0.size(); ++i)
      {
        if (min > count1s(M0[i]))
        {
          min = count1s(M0[i]);
          ind = i;
        }
      }

      // minimal solution
      for (size_t i = 0; i < primes.size(); ++i)
      {
        if (M0[ind] & (1ull << i))
        {
          int bit = 1 << var_cnt();
          int lit = 0;
          while (bit >>= 1)
          {
            if (!(primes[i].mask & bit))
            {
              used_var_flags[var_cnt() - lit - 1] = true;
            }
            ++lit;
          }
        }
      }

      return used_var_flags;
    }

    TYPE getTruthValue(const std::string &string)
    {
      std::string suffix_string;
      size_t index = string.find_first_of('=');
      if (index == std::string::npos)
        return 0;

      TYPE truth_value = 0;
      infix2Suffix(string.substr(index + 1), suffix_string);
      if (suffix_string == "1")
      {
        return ~truth_value;
      }
      else if (suffix_string == "0")
      {
        return truth_value;
      }

      std::vector<char> elements;  // index to vars or const or operator
      std::vector<char> operators; // 0: vars, 1:const, 2:operator
      std::string tmp_name;
      suffix_string += ' ';
      for (size_t i = 0; i < suffix_string.length(); i++)
      {
        char c = suffix_string.at(i);
        if (c == '+' || c == '*' || c == '~' || c == '@')
        {
          operators.push_back(2);
          elements.push_back(c);
        }
        else if (c == ' ' && tmp_name.length() > 0)
        {
          bool find = false;
          int j = 0;
          for (j = 0; j < var_cnt(); j++)
          {
            if (vars_[j] == tmp_name)
            {
              find = true;
              break;
            }
          }
          if (!find)
          {
            if (tmp_name == "0")
            {
              operators.push_back(1);
              elements.push_back(0);
            }
            else if (tmp_name == "1")
            {
              operators.push_back(1);
              elements.push_back(1);
            }
            else
            {
              assert(0);
              return 0;
            }
          }
          else
          {
            operators.push_back(0);
            elements.push_back(static_cast<char>(j));
          }
          tmp_name.clear();
        }
        else
        {
          if (c != ' ')
            tmp_name += c;
        }
      }

      std::vector<unsigned int> values(var_cnt()); // 0 or 1
      const int value_cnt = 1 << var_cnt();
      for (int i = 0; i < value_cnt; i++)
      {
        for (size_t j = 0; j < values.size(); j++)
          values[j] = 0;
        int bit = i;
        int index2 = 0;
        while (bit)
        {
          if (1 & bit)
          {
            values[index2] = 1;
          }
          else
          {
            values[index2] = 0;
          }
          bit = bit >> 1;
          index2++;
        }
        std::stack<unsigned int> stack;
        for (size_t j = 0; j < elements.size(); j++)
        {
          if (operators[j] == 0)
          { // vars
            stack.push(values[elements[j]]);
          }
          else if (operators[j] == 1)
          { // const 0, 1
            stack.push(elements[j]);
          }
          else
          {
            if (elements[j] == '~')
            {
              unsigned int v = stack.top();
              stack.pop();
              stack.push(static_cast<unsigned int>(1 - v));
            }
            else
            {
              // assert
              unsigned int v1 = stack.top();
              stack.pop();
              unsigned int v2 = stack.top();
              stack.pop();
              if (elements[j] == '+')
              {
                stack.push(v1 | v2);
              }
              else if (elements[j] == '*')
              {
                stack.push(v1 & v2);
              }
              else if (elements[j] == '@')
              {
                stack.push(v1 ^ v2);
              }
              else
              {
                assert(0);
              }
            }
          }
        }
        TYPE value = stack.top();
        truth_value |= (value << i);
      }
      return truth_value;
    }

  protected:
    int count1s(uint64_t x)
    {
      int cnt = 0;
      while (x)
      {
        cnt += static_cast<int>(x % 2);
        x >>= 1;
      }
      return cnt;
    }
    void mul(std::vector<uint64_t> &a, const std::vector<uint64_t> &b)
    {
      std::vector<uint64_t> v;
      for (size_t i = 0; i < a.size(); ++i)
      {
        for (size_t j = 0; j < b.size(); ++j)
        {
          v.push_back(a[i] | b[j]);
        }
      }
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      for (size_t i = 0; i < v.size() - 1; ++i)
      {
        for (size_t j = v.size() - 1; j > i; --j)
        {
          uint64_t z = v[i] & v[j];
          if ((z & v[i]) == v[i])
            v.erase(v.begin() + j);
          else if ((z & v[j]) == v[j])
          {
            uint64_t t = v[i];
            v[i] = v[j];
            v[j] = t;
            v.erase(v.begin() + j);
            j = v.size();
          }
        }
      }
      a = v;
    }

    std::string getName(Implicant &im)
    {
      std::string name;
      int bit = 1 << var_cnt();
      int lit = 0;
      bool f = false;
      while (bit >>= 1)
      {
        if (!(im.mask & bit))
        {
          if (f)
            name += "*";
          if (!(im.implicant & bit))
            name += "~";
          name += vars_[var_cnt() - lit - 1];
          f = true;
        }
        ++lit;
      }
      return name;
    }

    void infix2Suffix(const std::string &infix, std::string &suffix)
    {
      suffix.clear();
      std::stack<char> stack; // operator
      for (size_t i = 0; i < infix.length(); i++)
      {
        char c = infix.at(i);
        char temp;
        switch (c)
        {
        case ' ':
          break;
        case '(':
          stack.push(c);
          break;
        case '+':
          while (stack.size() != 0)
          {
            temp = stack.top();
            stack.pop();
            if (temp == '(')
            {
              stack.push(temp);
              break;
            }
            suffix += ' ';
            suffix += temp;
          }
          stack.push(c);
          suffix += ' ';
          break;
        case '@':
          while (stack.size() != 0)
          {
            temp = stack.top();
            stack.pop();
            if (temp == '(' || temp == '+')
            {
              stack.push(temp);
              break;
            }
            suffix += ' ';
            suffix += temp;
          }
          stack.push(c);
          suffix += ' ';
          break;
        case '*':
          while (stack.size() != 0)
          {
            temp = stack.top();
            stack.pop();
            if (temp == '(' || temp == '+' || temp == '@')
            {
              stack.push(temp);
              break;
            }
            suffix += ' ';
            suffix += temp;
          }
          stack.push(c);
          suffix += ' ';
          break;
        case '~':
          while (stack.size() != 0)
          {
            temp = stack.top();
            stack.pop();
            if (temp == '(' || temp == '+' || temp == '@' || temp == '*')
            {
              stack.push(temp);
              break;
            }
            suffix += ' ';
            suffix += temp;
          }
          stack.push(c);
          suffix += ' ';
          break;
        case ')':
          while (stack.size() != 0)
          {
            temp = stack.top();
            stack.pop();
            if (temp == '(')
            {
              break;
            }
            suffix += ' ';
            suffix += temp;
          }
          suffix += ' ';
          break;
        default:
          suffix += c;
          break;
        }
      }
      while (stack.size() != 0)
      {
        char temp = stack.top();
        suffix += ' ';
        suffix += temp;
        stack.pop();
      }
      return;
    }

  private:
    std::vector<std::string> vars_;
    std::string output_;
  };

  //====================================================================//
  //                           class pruner                             //
  //====================================================================//
  /**
   * @brief
   *
   * @tparam cell the cell will be ordered and pruned
   * @tparam compare the comparator of cell, not function pointer
   */
  template <typename Cell, typename Compare, typename Sizer, typename Valuer>
  class pruner
  {
  private:
    /**
     * @brief the default max store number for each cell set
     */
    int __store_num;

    /**
     *  the max cell number stored in each set
     */
    std::vector<short> __max_store_number;

    /**
     *  the ordered cells, index by cell size. In each list the cell will be unique, which means list will
     *  store only one cell if two or more cells are equivalent.
     */
    std::vector<std::set<Cell, Compare>> __cells;

    /**
     * the object of comparator
     */
    Compare __comparator;

    /**
     * @brief the object of sizer
     */
    Sizer __size;

    /**
     * @brief the object of __valuer
     */
    Valuer __valuer;

    /**
     * @brief the minimal value of pushed cells
     */
    float __min_value;

  public:
    /**
     * @brief Construct a new prune object
     *
     * @param max_cell_size the max size of cells
     * @param max_store_num the max number of cells stored for each size
     */
    pruner(int max_cell_size, int max_store_num = 4)
        : __store_num(max_store_num),
          __max_store_number(max_cell_size + 1, max_store_num),
          __cells(max_cell_size + 1),
          __min_value(100000000.0)
    {
    }

    ~pruner() {}

    /**
     * @brief push a cell into prune
     *
     * @param cell
     */
    bool push(const Cell &cell)
    {
      __min_value = __valuer(cell) < __min_value ? __valuer(cell) : __min_value;
      auto &set = __cells[__size(cell)];
      if (set.find(cell) != set.end())
        return false;
      set.emplace(cell);
      if (set.size() > __max_store_number[__size(cell)])
        set.erase(--set.end());
      return true;
    }

    /**
     * @brief collect the cells into given vector
     */
    void collect(std::vector<Cell> &cells, float comp_epsilon)
    {
      for (size_t i = 1; i != __cells.size(); ++i)
      {
        auto &cell_set = __cells[i];
        for (auto it = cell_set.rbegin(); it != cell_set.rend(); ++it)
        {
          // if ( cells.size() == 0 || __comparator( *it, cells.back(), loose ))
          if (cells.size() == 0 || __valuer(*it) < __valuer(cells.back()) + comp_epsilon)
          {
            cells.emplace_back(*it);
          }
        }
      }
      std::reverse(cells.begin(), cells.end());
    }

    /**
     * @brief collect the cells into given vector
     */
    void collect(std::vector<Cell> &cells, bool area_oriented)
    {
      // the max difference between min-cost cell and max-cost cell
      float diff = area_oriented ? 1.0 : 10.0;
      float value_upper = __min_value + diff;

      for (size_t i = 1; i != __cells.size(); ++i)
      {
        auto &cell_set = __cells[i];
        for (auto it = cell_set.rbegin(); it != cell_set.rend(); ++it)
        {
          if (area_oriented == false) // delay mode won't using area-size to prune cuts
          {
            cells.emplace_back(*it);
          }
          else if (__valuer(*it) <= value_upper && (cells.size() == 0 || __valuer(*it) < __valuer(cells.back())))
          {
            cells.emplace_back(*it);
          }
        }
      }
      std::reverse(cells.begin(), cells.end());
    }

    /**
     * @brief reset the max storing number of each set according to given vector
     *
     * @param new_max_store_num
     */
    void reset_store_num_upper(const std::vector<short> &new_max_store_num)
    {
      if (new_max_store_num.size() != __max_store_number.size())
      {
        std::cout << "reset prune size_upper failed.\n";
        return;
      }
      __max_store_number = new_max_store_num;
    }

    /**
     * @brief reset the max cell size to new_cell_size_upper and set the max store number by max_store_num(default = 4)
     *
     * @param new_cell_size_upper
     */
    void reset_cell_size_upper(int new_cell_size_upper, int new_store_num = 0)
    {
      __cells.clear();
      __cells.resize(new_cell_size_upper + 1);
      __max_store_number.clear();
      __store_num = new_store_num == 0 ? __store_num : new_store_num;
      __max_store_number.resize(new_cell_size_upper + 1, __store_num);
    }
  };

  //====================================================================//
  //                        class ordered_merge                         //
  //====================================================================//
  /**
   * @brief merging two ordered vector, no max merged size upper limit
   *
   * @tparam Value
   * @tparam _less the comparator " < " of Value
   * @param vec_a
   * @param vec_b
   * @param result
   * @return size_t
   */
  template <class Value, class _less>
  void ordered_merge(const std::vector<Value> &vec_a,
                     const std::vector<Value> &vec_b,
                     std::vector<Value> &result)
  {
    size_t idx_a = 0;
    size_t idx_b = 0;
    _less less;

    result.reserve(vec_a.size() + vec_a.size());

    while (idx_a != vec_a.size() && idx_b != vec_b.size())
    {
      if (less(vec_a[idx_a], vec_b[idx_b]))
      {
        result.emplace_back(vec_b[idx_b++]);
      }
      else if (less(vec_b[idx_b], vec_a[idx_a]))
      {
        result.emplace_back(vec_a[idx_a++]);
      }
      else
      {
        result.emplace_back(vec_a[idx_a++]);
        ++idx_b;
      }
    }

    if (idx_a == vec_a.size())
    {
      while (idx_b != vec_b.size())
      {
        result.emplace_back(vec_b[idx_b++]);
      }
      return;
    }

    if (idx_b == vec_b.size())
    {
      while (idx_a != vec_a.size())
      {
        result.emplace_back(vec_a[idx_a++]);
      }
    }
  }

  /**
   * @brief ordered vector merge with size uppper limit
   *
   * @tparam Value
   * @tparam comp_value the property of Value used to compare, must overload < !=
   * @param vec_a
   * @param vec_b
   * @param result
   * @param size_upper
   * @return true the merged size has greater than size upper
   * @return false
   */
  template <class Value, class comp_value>
  bool ordered_merge(const std::vector<Value> &vec_a,
                     const std::vector<Value> &vec_b,
                     std::vector<Value> &result,
                     size_t size_upper)
  {
    comp_value val;

    // result.resize(size_upper);
    result.reserve(size_upper);

    size_t sizeA = vec_a.size();
    size_t sizeB = vec_b.size();

    // both vector are the largest
    if (sizeA == size_upper && sizeB == size_upper)
    {
      for (int i = 0; i != sizeA; ++i)
      {
        if (val(vec_a[i]) != val(vec_b[i]))
        {
          return true;
        }
        result.emplace_back(vec_a[i]);
        // result[i] = vec_a[i];
      }
      return false;
    }

    int i, k, c;
    i = k = c = 0;

    while (true)
    {
      if (val(vec_a[i]) < val(vec_b[k]))
      {
        // result[c++] = vec_a[i++];
        result.emplace_back(vec_a[i++]);
        if (i == sizeA)
        {
          goto FlushVec1;
        }
      }
      else if (val(vec_b[k]) < val(vec_a[i]))
      {
        // result[c++] = vec_b[k++];
        result.emplace_back(vec_b[k++]);
        if (k == sizeB)
        {
          goto FlushVec0;
        }
      }
      else
      {
        // result[c++] = vec_a[i++];
        result.emplace_back(vec_a[i++]);
        ++k;
        if (i == sizeA)
        {
          goto FlushVec1;
        }
        if (k == sizeB)
        {
          goto FlushVec0;
        }
      }
    }

  FlushVec0:
    c = result.size();
    if (c + sizeA > size_upper + i)
    {
      return true;
    }
    while (i < sizeA)
    {
      result.emplace_back(vec_a[i++]);
      // result[c++] = vec_a[i++];
    }
    return false;

  FlushVec1:
    c = result.size();
    if (c + sizeB > size_upper + k)
    {
      return true;
    }
    while (k < sizeB)
    {
      result.emplace_back(vec_b[k++]);
      // result[c++] = vec_b[k++];
    }
    return false;
  }

  //====================================================================//
  //                           truth table                              //
  //====================================================================//
  class truth_table
  {
  public:
    truth_table() {}
    truth_table(const std::string &expression) : expression_(expression) {}
    truth_table(const truth_table &tt) = default;

    /**
     * @brief truth table AND
     *
     * @param tt
     * @return truth_table result new truth table
     */
    truth_table operator&(const truth_table &tt)
    {
      std::string exp = this->expression_;
      exp = "(" + exp + ")*(" + tt.expression_ + ")";
      return truth_table(exp);
    }

    /**
     * @brief truth table AND
     *
     * @param tt
     * @return truth_table result new truth table
     */
    truth_table &operator&=(const truth_table &tt)
    {
      // std::string exp = this->expression_;
      // exp = "(" + exp + ")&(" + tt.expression_ + ")";
      // return truth_table(exp);
      if (expression_.size() == 0)
        expression_ = tt.expression_;
      else
      {
        expression_ = "(" + expression_ + ")*(" + tt.expression_ + ")";
      }
      return *this;
    }

    /**
     * @brief using argument to init tt
     *
     * @param expresstion  argument name
     * @return truth_table&
     */
    truth_table &operator=(const std::string &expresstion)
    {
      expression_ = expresstion;
      return *this;
    }

    /**
     * @brief truth table OR
     *
     * @param tt
     * @return truth_table result new truth table
     */
    truth_table operator|(const truth_table &tt)
    {
      std::string exp = this->expression_;
      exp = "(" + exp + ")|(" + tt.expression_ + ")";
      return truth_table(exp);
    }

    /**
     * @brief truth table INV
     *
     * @return truth_table result new truth table
     */
    truth_table operator~()
    {
      std::string exp = this->expression_;
      exp = "~(" + exp + ")";
      return truth_table(exp);
    }

    /**
     * @brief return the funtion by sop expression type
     *
     * @return const std::string& expression
     */
    const std::string &function() const { return expression_; }

    /**
     * @brief clear the truth table
     */
    void clear() { expression_.clear(); }

    /**
     * @brief swap two truth table
     */
    void swap(truth_table &tt) { expression_.swap(tt.expression_); }

  private:
    std::string expression_;
  };

  //===----------------------------------------------------------------------===//
  //                               usful functions                              //
  //===----------------------------------------------------------------------===//
  static std::string int2Hex(unsigned long val, int len)
  {
    std::string init;
    init.resize(len);
    for (int i = 0; i < len; ++i)
    {
      int cc = val & 15UL;
      if (cc < 10)
        cc += static_cast<int>('0');
      else
        cc += static_cast<int>('A') - 10;
      init[len - i - 1] = static_cast<char>(cc);
      val >>= 4;
    }
    return init;
  }

} // end namespace