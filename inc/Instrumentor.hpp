//
// Basic instrumentation profiler by Cherno

// Usage: include this header file somewhere in your code (eg. precompiled
// header), and then use like:
//
// Instrumentor::Get().BeginSession("Session Name");        // Begin session
// {
//     InstrumentationTimer timer("Profiled Scope Name");   // Place code like
//     this in scopes you'd like to include in profiling
//     // Code
// }
// Instrumentor::Get().EndSession();                        // End Session
//
// You will probably want to macro-fy this, to switch on/off easily and use
// things like __FUNCSIG__ for the profile name.
//
#pragma once

#include <algorithm>
#include <chrono>
#include <fstream>
#include <string>
#include <thread>

#include "mpi.h"

struct ProfileResult {
  std::string Name;
  long long Start, End;
  uint32_t ThreadID;
  int ProcID;
};

struct InstrumentationSession {
  std::string Name;
};

class Instrumentor {
 private:
  InstrumentationSession *m_CurrentSession;
  std::ofstream m_OutputStream;
  int m_ProfileCount;

 public:
  Instrumentor() : m_CurrentSession(nullptr), m_ProfileCount(0) {}

  void BeginSession(const std::string &
                        name)  //, const std::string& filepath = "results.json")
  {
    m_OutputStream.open(name);
    WriteHeader();
    m_CurrentSession = new InstrumentationSession{name};
  }

  void EndSession() {
    WriteFooter();
    m_OutputStream.close();
    delete m_CurrentSession;
    m_CurrentSession = nullptr;
    m_ProfileCount = 0;
  }

  void WriteProfile(const ProfileResult &result) {
    if (m_ProfileCount++ > 0) m_OutputStream << ",";

    std::string name = result.Name;
    std::replace(name.begin(), name.end(), '"', '\'');

    m_OutputStream << "{";
    m_OutputStream << "\"cat\":\"function\",";
    m_OutputStream << "\"dur\":" << (result.End - result.Start) << ',';
    m_OutputStream << "\"name\":\"" << name << "\",";
    m_OutputStream << "\"ph\":\"X\",";
    m_OutputStream << "\"pid\":" << result.ProcID << ",";
    m_OutputStream << "\"tid\":" << result.ThreadID << ",";
    m_OutputStream << "\"ts\":" << result.Start;
    m_OutputStream << "}";

    m_OutputStream.flush();
  }

  void WriteHeader() {
    m_OutputStream << "{\"otherData\": {},\"traceEvents\":[";
    m_OutputStream.flush();
  }

  void WriteFooter() {
    m_OutputStream << "]}";
    m_OutputStream.flush();
  }

  static Instrumentor &Get() {
    static Instrumentor instance;
    return instance;
  }
};

class InstrumentationTimer {
 public:
  explicit InstrumentationTimer(const char *name)
      : m_Name(name),
        m_StartTimepoint(std::chrono::high_resolution_clock::now()),
        m_Stopped(false) {}

  ~InstrumentationTimer() {
    if (!m_Stopped) Stop();
  }

  void Stop() {
    auto endTimepoint = std::chrono::high_resolution_clock::now();

    long long start = std::chrono::time_point_cast<std::chrono::microseconds>(
                          m_StartTimepoint)
                          .time_since_epoch()
                          .count();
    long long end =
        std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint)
            .time_since_epoch()
            .count();

    uint32_t threadID =
        std::hash<std::thread::id>{}(std::this_thread::get_id());
    int myProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
    Instrumentor::Get().WriteProfile({m_Name, start, end, threadID, myProc});

    m_Stopped = true;
  }

 private:
  const char *m_Name;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
  bool m_Stopped;
};

#define PROFILING 0
#if PROFILING
#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)
#define PROFILE_FUNCTION(name) PROFILE_SCOPE(__FUNCTION__)
#else
#define PROFILE_SCOPE(name)
#define PROFILE_FUNCTION(name) PROFILE_SCOPE(__FUNCTION__)
#endif
