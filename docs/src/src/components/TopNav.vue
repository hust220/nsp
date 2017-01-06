<template>
  <el-row class="top-nav">
    <el-menu default-active="1" class="el-menu-demo" mode="horizontal" @select="handleSelect">
      <el-col :span="18" :offset="3">
        <el-menu-item index="1">
          <router-link to="/home">Home</router-link>
        </el-menu-item>
        <el-menu-item index="1">
          <router-link to="/doc">Documentation</router-link>
        </el-menu-item>
        <el-menu-item index="1">
          <router-link to="/api">API</router-link>
        </el-menu-item>
        <el-menu-item index="1">
          <router-link to="/download">Download</router-link>
        </el-menu-item>
        <el-menu-item index="1">
        <el-menu-item index="3"><a href="mailto:wj_hust08@hust.edu.cn?cc=yxiao@hust.edu.cn">Contact</a></el-menu-item>
      </el-col>
    </el-menu>
    <div class="line"></div>
  </el-row>
</template>

<script>
  export default {
    name: 'nav',
    methods: {
      handleSelect(key, keyPath) {
        // console.log(key, keyPath)
      }
    }
  }
</script>

<style scoped>
a {
  color: rgb(71, 86, 105);
  text-decoration: none;
}
</style>