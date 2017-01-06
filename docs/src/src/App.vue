<template>
  <div id="app">
    <!-- Top -->
    <el-row>
      <el-col :span="18" :offset="3">
        <top></top>
      </el-col>
    </el-row>

    <!-- TopNav-->
    <top-nav></top-nav>

    <!-- Page -->
    <el-row>
      <el-col :span="18" :offset="3">
        <router-view></router-view>
      </el-col>
    </el-row>

  </div>
</template>

<script>
  import Top from './components/Top'
  import TopNav from './components/TopNav'
  import { bus } from './bus.js'

  export default {
    name: 'app',
    components: {
      'top-nav': TopNav,
      'top': Top,
    }
  }
</script>

<style>
#app {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
}
body {
  margin: 0px;
}
</style>